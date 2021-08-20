import h5py
import logging
import operator
import numpy as np

import westpa
from westpa.cli.core import w_init
from westpa.cli.core import w_run
from westpa.core.extloader import get_object
from westpa.core.segment import Segment

import json

import os
import shutil
import pickle
import importlib.util

import tqdm

import mdtraj as md
from rich.logging import RichHandler

# Ensure this is installed via pip. msm_we's setup.py is all set up for that.
# Navigate to the folder where msm_we is, and run python3 -m pip install .
#   If you're doing development on msm_we, add the -e flag to pip, i.e. "python3 -m pip install -e ."
#   -e will install it in editable mode, so changes to msm_we will take effect next time it's imported.
#   Otherwise, if you modify the msm_we code, you'll need to re-install it through pip.
from msm_we import msm_we

EPS = np.finfo(np.float64).eps

log = logging.getLogger(__name__)
log.setLevel("INFO")
log.propagate = False
log.addHandler(RichHandler())

msm_we_logger = logging.getLogger("msm_we.msm_we")
msm_we_logger.setLevel("INFO")

# Map structure types to extensions.
# This tells the plugin what extension to put on generated start-state files.
STRUCT_EXTENSIONS = {
    md.formats.PDBTrajectoryFile: "pdb",
    md.formats.AmberRestartFile: "rst7",
}


def check_target_reached(h5_filename):
    """
    Check if the target state was reached, given the data in a WEST H5 file.

    Parameters
    ----------
    h5_filename: string
        Path to a WESTPA HDF5 data file
    """
    with h5py.File(h5_filename, 'r') as h5_file:
        # Get the key to the final iteration. Need to do -2 instead of -1 because there's an empty-ish final iteration
        #   written.
        for iteration_key in list(h5_file['iterations'].keys())[-2:0:-1]:
            endpoint_types = h5_file[f'iterations/{iteration_key}/seg_index']['endpoint_type']
            if Segment.SEG_ENDPOINT_RECYCLED in endpoint_types:
                return True
    return False


def fix_deprecated_initialization(initialization_state):
    """
    I changed my initialization JSON schema to use underscores instead of hyphens so I can directly expand it into
    keywords arguments to w_init. This just handles any old-style JSON files I still had, so they don't choke and die.
    """

    log.info(f"Starting processing, dict is now {initialization_state}")

    # Some of my initial files had this old-style formatting. Handle it for now, but remove eventually
    for old_key, new_key in [
        ('tstate-file', 'tstate_file'),
        ('bstate-file', 'bstate_file'),
        ('sstate-file', 'sstate_file'),
        ('segs-per-state', 'segs_per_state'),
    ]:
        if old_key in initialization_state.keys():
            log.warning(
                f"This initialization JSON file uses the deprecated  " f"hyphenated form for  {old_key}. Replace with underscores."
            )

            value = initialization_state.pop(old_key)
            initialization_state[new_key] = value

    log.info(f"Finished processing, dict is now {initialization_state}")
    return initialization_state


# TODO: Break this out into a separate module, let it be specified (if it's necessary) as a plugin option
#   This may not always be required -- i.e. you may be able to directly output to the h5 file in your propagator
def prepare_coordinates(plugin_config, h5file, we_h5filename):
    """
    Copy relevant coordinates from trajectory files into <iteration>/auxdata/coord of the h5 file.

    Directly modifies the input h5 file.

    Adds ALL coordinates to auxdata/coord.

    Adapted from original msmWE collectCoordinates.py script.

    Parameters
    ----------
    plugin_config: YAMLConfig object
        Stores the configuration options provided to the plugin in the WESTPA configuration file

    h5file: h5py.File
        WESTPA h5 data file

    we_h5filename: string
        Name of the WESTPA h5 file
    """

    refPDBfile = plugin_config.get('ref_pdb_file')
    modelName = plugin_config.get('model_name')

    # TODO: Don't need this explicit option, use WEST_SIM_ROOT or something
    WEfolder = plugin_config.get('we_folder')

    parentTraj = plugin_config.get('parent_traj_filename')
    childTraj = plugin_config.get('child_traj_filename')

    model = msm_we.modelWE()
    log.info('Preparing coordinates...')

    # Only need the model to get the number of iterations and atoms
    # TODO: Replace this with something more lightweight, get directly from WE
    model.initialize(we_h5filename, refPDBfile, modelName)
    model.get_iterations()

    n_iter = None
    for n_iter in tqdm.tqdm(range(1, model.maxIter)):

        nS = model.numSegments[n_iter - 1].astype(int)
        coords = np.zeros((nS, 2, model.nAtoms, 3))
        dsetName = "/iterations/iter_%08d/auxdata/coord" % int(n_iter)

        coords_exist = False
        try:
            dset = h5file.create_dataset(dsetName, np.shape(coords))
        except (RuntimeError, ValueError):
            log.debug('coords exist for iteration ' + str(n_iter) + ' NOT overwritten')
            coords_exist = True
            continue

        for iS in range(nS):
            trajpath = WEfolder + "/traj_segs/%06d/%06d" % (n_iter, iS)

            try:
                coord0 = np.squeeze(md.load(f'{trajpath}/{parentTraj}', top=model.reference_structure.topology)._xyz)
            except OSError:
                log.warning("Parent traj file doesn't exist, loading reference structure coords")
                coord0 = np.squeeze(model.reference_structure._xyz)

            coord1 = np.squeeze(md.load(f'{trajpath}/{childTraj}', top=model.reference_structure.topology)._xyz)

            coords[iS, 0, :, :] = coord0
            coords[iS, 1, :, :] = coord1

        if not coords_exist:
            dset[:] = coords

    log.debug(f"Wrote coords for {n_iter} iterations.")


def msmwe_compute_ss(plugin_config, west_files, last_iter):
    """
    Prepare and initialize an msm_we model, and use it to predict a steady-state distribution.

    1. Load coordinate data
    2. Perform dimensionality reduction
    3. Compute flux and transition matrices
    4. Compute steady-state distribution (via eigenvectors of transition matrix)
    5. Compute target-state flux

    TODO
    ----
    This function does far too many things. Break it up a bit.

    Parameters
    ----------
    plugin_config: YAMLConfig object
        Stores the configuration options provided to the plugin in the WESTPA configuration file

    last_iter: int
        The last WE iteration to use for computing steady-state.

    Returns
    -------
    ss_alg: np.ndarray
        The steady-state distribution

    ss_flux: float
        Flux into target state

    model: modelWE object
        The modelWE object produced for analysis.
    """

    n_lag = 0

    log.debug("Initializing msm_we")

    # TODO: Refactor this to use westpa.core.extloader.get_object
    #   I'm reinventing the wheel a bit here, I can replace almost all this code w/ that
    # ##### Monkey-patch modelWE with the user-override functions
    override_file = plugin_config.get('user_functions')

    # First, import the file with the user-override functions
    # This is a decently janky implementation, but it seems to work, and I don't know of a better way of doing it.
    # This is nice because it avoids mucking around with the path, which I think is a Good Thing.

    # We're given a path to the user-specified file containing overrides
    # This comes from https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path

    # I don't think the name provided here actually matters
    user_override_spec = importlib.util.spec_from_file_location("override_module", override_file)
    user_overrides = importlib.util.module_from_spec(user_override_spec)

    # Make the functions that were overriden in override_file available in the namespace under user_overrides
    user_override_spec.loader.exec_module(user_overrides)

    # So now we can do the actual monkey-patching of modelWE.
    # We monkey-patch at the module level rather than just override the function in the instanced object
    #   so that the functions retain access to self.
    msm_we.modelWE.processCoordinates = user_overrides.processCoordinates
    # ##### Done with monkey-patching.

    model = msm_we.modelWE()

    fileSpecifier = ' '.join(west_files)
    refPDBfile = plugin_config.get('ref_pdb_file')
    modelName = plugin_config.get('model_name')
    n_clusters = plugin_config.get('n_clusters')
    tau = plugin_config.get('tau', None)

    # Fire up the model object
    # (Eventually this will just go in __init__)

    # In RestartXX/RunYY fileSpecifier is a list of all Restart{0..XX}/Run{1..YY}/west.h5
    # model.initialize(fileSpecifier, refPDBfile, initPDBfile, modelName)
    model.initialize(fileSpecifier, refPDBfile, modelName, tau)

    # Set some model parameters
    WEtargetp1 = plugin_config.get('target_pcoord1', None)
    if WEtargetp1 is not None:
        log.warning(
            "Using a deprecated target pcoord specification. "
            "You should change this to explicitly give target_pcoord1_min and target_pcoord1_max as bounds. "
            "Defaulting min to -np.inf for now.."
        )

    model.WEtargetp1_bounds = [
        plugin_config.get('target_pcoord1_min', -np.inf),
        plugin_config.get('target_pcoord1_max', WEtargetp1),
    ]
    model.WEbasisp1_bounds = [plugin_config.get('basis_pcoord1_min'), plugin_config.get('basis_pcoord1_max')]

    model.pcoord_ndim0 = plugin_config.get('pcoord_ndim0')
    model.dimReduceMethod = plugin_config.get('dim_reduce_method')

    model.n_lag = n_lag

    # last_iter_cluster = last_iter  # model.maxIter-1 #last iter often not complete
    # i = last_iter_cluster
    # coordSet = np.zeros((0, model.nAtoms, 3))  # extract coordinate libraries for clustering

    log.debug("Loading in iteration data.. (this could take a while)")
    # log.debug(f'coord shape is {coordSet.shape}')

    # First dimension is the total number of segments
    model.get_iterations()

    model.get_coordSet(last_iter)

    model.dimReduce()

    first_iter, last_iter = model.first_iter, model.last_iter

    clusterFile = modelName + "_clusters_s" + str(first_iter) + "_e" + str(last_iter) + "_nC" + str(n_clusters) + ".h5"
    # TODO: Uncomment this to actually load the clusterFile if it exists.  For now, disable for development.
    exists = os.path.isfile(clusterFile)
    exists = False
    log.warning("Skipping any potential cluster reloading!")

    # If a cluster file with the name corresponding to these parameters exists, load clusters from it.
    if exists:
        log.debug("loading clusters...")
        model.load_clusters(clusterFile)
    # Otherwise, do the clustering (which will create and save to that file)
    else:
        # FIXME: This gives the wrong shape, but loading from the clusterfile gives the right shape
        log.debug("clustering coordinates into " + str(n_clusters) + " clusters...")
        model.cluster_coordinates(n_clusters)

    first_iter = 1
    model.get_fluxMatrix(n_lag, first_iter, last_iter)  # extracts flux matrix, output model.fluxMatrixRaw
    log.debug(f"Unprocessed flux matrix has shape {model.fluxMatrixRaw.shape}")
    model.organize_fluxMatrix()  # gets rid of bins with no connectivity, sorts along p1, output model.fluxMatrix
    model.get_Tmatrix()  # normalizes fluxMatrix to transition matrix, output model.Tmatrix

    log.debug(f"Processed flux matrix has shape {model.fluxMatrix.shape}")

    model.get_steady_state_algebraic()  # gets steady-state from eigen decomp, output model.pSS
    model.get_steady_state_target_flux()  # gets steady-state target flux, output model.JtargetSS

    # Why is model.pss sometimes the wrong shape? It's "occasionally" returned as a nested array.
    # Squeeze fixes it and removes the dimension of length 1, but why does it happen in the first place?

    if type(model.pSS) is np.matrix:
        ss_alg = np.squeeze(model.pSS.A)
    else:
        ss_alg = np.squeeze(model.pSS)

    ss_flux = model.JtargetSS

    log.debug("Got steady state:")
    log.debug(ss_alg)
    log.debug(ss_flux)

    log.info("Completed flux matrix calculation and steady-state estimation!")

    return ss_alg, ss_flux, model


class RestartDriver:
    """
    WESTPA plugin to automatically handle estimating steady-state from a WE run, re-initializing a new WE run in that
    steady-state, and then running that initialized WE run.

    Data from the previous run will be stored in the restart<restart_number>/ subdirectory of $WEST_SIM_ROOT.

    This plugin depends on having the start-states implementation in the main WESTPA code, which allows initializing
    a WE run using states that are NOT later used for recycling.

    These are used so that when the new WE run is initialized, initial structure selection is chosen by w_init, using
    weights assigned to the start-states based on MSM bin weight and WE segment weight.

    Since it closes out the current WE run and starts a new one, this plugin should run LAST, after all other plugins.
    """

    def __init__(self, sim_manager, plugin_config):
        """
        Initialize the RestartDriver plugin.

        Pulls the data_manager and sim_manager from the WESTPA run that just completed, along with
        """

        westpa.rc.pstatus("Restart plugin initialized")

        if not sim_manager.work_manager.is_master:
            westpa.rc.pstatus("Reweighting not master, skipping")
            return

        self.data_manager = sim_manager.data_manager
        self.sim_manager = sim_manager

        self.plugin_config = plugin_config

        self.restart_file = plugin_config.get('restart_file', 'restart.dat')
        self.initialization_file = plugin_config.get('initialization_file', 'restart_initialization.json')

        self.coord_len = plugin_config.get('coord_len', 2)
        self.n_restarts = plugin_config.get('n_restarts', -1)
        self.n_runs = plugin_config.get('n_runs', 1)

        # .get() might return this as a bool anyways, but be safe
        self.debug = bool(plugin_config.get('debug', False))
        if self.debug:
            log.setLevel("DEBUG")
            msm_we_logger.setLevel("DEBUG")

        # Default to using all restarts
        self.restarts_to_use = plugin_config.get('n_restarts_to_use', self.n_restarts)
        assert self.restarts_to_use > 0 or self.restarts_to_use == -1, "Invalid number of restarts to use"
        if self.restarts_to_use >= 1:
            assert (
                self.restarts_to_use == self.restarts_to_use // 1
            ), "If choosing a decimal restarts_to_use, must be between 0 and 1."

        struct_filetype = plugin_config.get('struct_filetype', 'mdtraj.formats.PDBTrajectoryFile')
        self.struct_filetype = get_object(struct_filetype)

        # This should be low priority, because it closes the H5 file and starts a new WE run. So it should run LAST
        #   after any other plugins.
        self.priority = plugin_config.get('priority', 100)  # I think a big number is lower priority...

        sim_manager.register_callback(sim_manager.finalize_run, self.prepare_new_we, self.priority)

        # Initialize data
        self.ss_alg = None
        self.ss_dist = None
        self.model = None

    def get_original_bins(self):
        """
        Obtains the WE bins and their probabilities at the end of the previous iteration.

        Returns
        -------
        bins : np.ndarray
            Array of WE bins

        binprobs: np.ndarray
            WE bin weights
        """

        we_driver = self.sim_manager.we_driver
        bins = we_driver.next_iter_binning
        n_bins = len(bins)
        binprobs = np.fromiter(map(operator.attrgetter('weight'), bins), dtype=np.float64, count=n_bins)

        return bins, binprobs

    @property
    def cur_iter(self):
        """
        Get the current WE iteration.

        Returns
        -------
        int: The current iteration. Subtract one, because in finalize_run the iter has been incremented
        """
        return self.sim_manager.n_iter - 1

    @property
    def is_last_iteration(self):
        """
        Get whether this is the last iteration in this WE run.

        Returns
        -------
        bool: Whether the current iteration is the final iteration
        """

        final_iter = self.sim_manager.max_total_iterations

        return self.cur_iter == final_iter

    def prepare_new_we(self):
        """
        This function prepares a new WESTPA simulation using haMSM analysis to accelerate convergence.

        The marathon functionality does re-implement some of the functionality of w_multi_west.
        However, w_multi_west merges independent WE simulations, which may or may not be desirable.
        I think for the purposes of this, it's good to keep the runs completely independent until haMSM model building.
        Either that, or I'm just justifying not having known about w_multi_west when I wrote this. TBD.

        # TODO: Replace all manual path-building with pathlib

        The algorithm is as follows:
            1. Check to see if we've just completed the final iteration
            2. Handle launching multiple runs, if desired
            2. Build haMSM
            3. Obtain structures for each haMSM bin
            4. Make each structure a start-state, with probability set by (MSM-bin SS prob / # structures in bin)
            5. Potentially some renormalization?
            6. Start new WE simulation
        """

        # Do nothing if it's not the final iteration
        if not self.is_last_iteration:
            print(self.cur_iter)
            return

        log.debug("Final iteration, preparing restart")

        restart_state = {'restarts_completed': 0, 'runs_completed': 0}

        # TODO: Check for the existence of the extension lockfile here
        doing_extension = False

        # Look for a restart.dat file to get the current state (how many restarts have been performed already)
        if os.path.exists(self.restart_file):
            with open(self.restart_file, 'r') as fp:
                restart_state = json.load(fp)

        # This is the final iteration of a run, so mark this run as completed
        restart_state['runs_completed'] += 1

        # Make the folder to store data for this marathon
        restart_directory = f"restart{restart_state['restarts_completed']}"
        run_directory = f"{restart_directory}/run{restart_state['runs_completed']}"
        if not os.path.exists(run_directory):
            os.makedirs(run_directory)

        # Write coordinates to h5
        prepare_coordinates(self.plugin_config, self.data_manager.we_h5file, self.data_manager.we_h5filename)

        for data_folder in ['traj_segs', 'seg_logs']:
            old_path = data_folder
            new_path = f"{run_directory}/{old_path}"

            log.debug(f"Moving {old_path} to {new_path}")

            if os.path.exists(new_path):
                log.error(f"{new_path} already exists. Removing and overwriting.")
                shutil.rmtree(new_path)

            try:
                os.rename(old_path, new_path)
            except FileNotFoundError:
                log.warning(f"Folder {old_path} was not found." "This may be normal, but check your configuration.")
            else:
                # Make a new data folder for the next run
                os.mkdir(old_path)

        last_run = restart_state['runs_completed'] >= self.n_runs
        last_restart = restart_state['restarts_completed'] >= self.n_restarts

        # We've just finished a run. Let's check if we have to do any more runs in this marathon before doing a restart.
        #   In the case of n_runs == 1, then we're just doing a single run and restarting it every so often.
        #   Otherwise, a marathon consists of multiple runs,  and restarts are performed between marathons.
        if last_run:
            log.info(f"All {self.n_runs} runs in this marathon completed.")

            if last_restart:
                log.info("All restarts completed! Performing final analysis.")

            else:
                log.info("Proceeding to prepare a restart.")

                # Duplicating this is  gross, but given the structure here, my options are either put it above these ifs
                #   entirely, meaning it'll be unnecessarily run at the end of the final restart, or duplicate it below.
                log.info("Preparing coordinates for this run.")

                # Now, continue on to haMSM calculation below.

        # If we have more runs left to do in this marathon, prepare them
        elif not last_run:

            log.info(f"Run {restart_state['runs_completed']}/{self.n_runs} completed.")

            # TODO: Initialize a new run, from the same configuration as this run was
            #   On the 1st run, I can write bstates/tstates/sstates into restart files, and use those for spawning
            #   subsequent runs in the marathon. That way, I don't make unnecessary copies of all those.
            # Basis and target states are unchanged. Can I get the original parameters passed to w_init?
            # Ideally, I should be able to call w_init with the exact same parameters that went to it the first time
            initialization_state = {
                'tstate_file': None,
                'bstate_file': None,
                'sstate_file': None,
                'tstates': None,
                'bstates': None,
                'sstates': None,
                'segs_per_state': None,
            }

            # TODO: Implement this, and get rid of the initialization_file usage right below. Placeholder for now.
            if restart_state['runs_completed'] == 1:

                # Get and write basis, target, start states and segs per state for this marathon to disk
                pass

            # Save the WESTPA h5 data from this run
            self.data_manager.finalize_run()
            shutil.copyfile('west.h5', f"{run_directory}/west.h5")

            if not doing_extension:
                if os.path.exists(self.initialization_file):
                    with open(self.initialization_file, 'r') as fp:
                        initialization_dict = json.load(fp)
                        initialization_dict = fix_deprecated_initialization(initialization_dict)
                        initialization_state.update(initialization_dict)
                else:
                    raise Exception(
                        "No initialization JSON file provided -- " "I don't know how to start new runs in this marathon."
                    )

                westpa.rc.pstatus(
                    f"\n\n===== Restart {restart_state['restarts_completed']}, "
                    + f"Run {restart_state['runs_completed']+1} initializing =====\n"
                )

                westpa.rc.pstatus(
                    f"\nRun: \n\t w_init --tstate-file {initialization_state['tstate_file']} "
                    + f"--bstate-file {initialization_state['bstate_file']} "
                    f"--sstate-file {initialization_state['sstate_file']} "
                    f"--segs-per-state {initialization_state['segs_per_state']}\n"
                )

                w_init.initialize(
                    **initialization_state, shotgun=False,
                )

            # If we're doing an extension set
            elif doing_extension:

                # TODO:  Instead of w_initting a new iteration, copy the files from restart0/runXX back into ./

                # Then, continue as normal
                pass

            with open(self.restart_file, 'w') as fp:
                json.dump(restart_state, fp)

            log.info("New WE run ready!")
            westpa.rc.pstatus(
                f"\n\n===== Restart {restart_state['restarts_completed']}, "
                + f"Run {restart_state['runs_completed']+1} running =====\n"
            )

            w_run.run_simulation()
            return

        log.debug(f"{restart_state['restarts_completed']}/{self.n_restarts} restarts completed")

        # Build the haMSM
        log.debug("Initializing haMSM")

        # Need to write the h5 file and close it out, but I need to get the current bstates first.
        original_bstates = self.sim_manager.current_iter_bstates
        if original_bstates is None:
            original_bstates = self.data_manager.get_basis_states(self.sim_manager.n_iter - 1)

        assert original_bstates is not None, "Bstates are none in the current iteration"

        original_tstates = self.data_manager.get_target_states(self.cur_iter)

        # Flush h5 file writes and copy it to the run directory
        self.data_manager.finalize_run()
        shutil.copyfile(self.data_manager.we_h5filename, f"{run_directory}/west.h5")

        log.debug("Building haMSM and computing steady-state")

        # Use all files in all restarts
        # Restarts index at 0, because there's  a 0th restart before you've... restarted anything.
        # Runs index at 1, because Run 1 is the first run.
        # TODO: Let the user pick last half or something in the plugin config.
        marathon_west_files = []
        # When doing the first restart, restarts_completed is 0 (because the first restart isn't complete yet) and
        #   the data generated during this restart is in /restart0.
        # So when doing the Nth restart, restarts_completed is N-1

        # If set to -1, use all restarts
        if self.restarts_to_use == -1:
            last_N_restarts = 1 + restart_state['restarts_completed']
        # If this is an integer, use the last N restarts
        elif self.restarts_to_use >= 1:
            last_N_restarts = self.restarts_to_use
        # If it's a decimal between 0 and 1, use it as a fraction
        # At restart 1, and a fraction of 0.5, this should just use restart 1
        elif 0 < self.restarts_to_use < 1:
            last_N_restarts = int(self.restarts_to_use * (1 + restart_state['restarts_completed']))

            # If this fraction is <1, use all until it's not
            if last_N_restarts < 1:
                last_N_restarts = 1 + restart_state['restarts_completed']

        log.debug(f"Last N is {last_N_restarts}")
        first_restart = max(1 + restart_state['restarts_completed'] - last_N_restarts, 0)
        usable_restarts = range(first_restart, 1 + restart_state['restarts_completed'])

        log.info(
            f"At restart {restart_state['restarts_completed']}, building haMSM using data from restarts {list(usable_restarts)}"
        )
        for restart_number in usable_restarts:
            for run_number in range(1, 1 + restart_state['runs_completed']):

                west_file_path = f"restart{restart_number}/run{run_number}/west.h5"
                marathon_west_files.append(west_file_path)

        log.debug(f"WESTPA datafile for analysis are {marathon_west_files}")

        #
        # If this is the first restart, check to see if you got any target state flux
        if restart_state['restarts_completed'] == 0:
            pass

            # TODO: Check to see if you got any target flux in ANY runs
            target_reached = True

            # If you reached the target, clean up from the extensions and then continue as normal
            if target_reached:

                # TODO: Remove the doing_extensions.lck lockfile

                # TODO: Restore the original, non-extended west.cfg

                # Then continue as normal
                pass

            # If no runs reached the target, then we need to extend them
            elif not target_reached:

                # TODO: Write the doing_extensions.lck "lockfile" to indicate we're in extend mode (or keep if exists)

                # TODO: Reset runs_completed to 0

                # TODO: Copy traj_segs, seg_logs, and west.h5 for restart0/run0 back into ./
                #       Later: (May only need to copy the latest iteration traj_segs, to avoid tons of back and forth)

                # TODO: Backup (copy) west.cfg to west_normal.cfg  (if west_normal.cfg does not exist)
                #       If it DOES exist, then we must be extending again, so extend on the current west.cfg

                # TODO: Create a new west.cfg by extending west.propagation.max_total_iterations

                # TODO: w_run (which will use the extended west.cfg)
                pass

        ss_dist, ss_flux, model = msmwe_compute_ss(self.plugin_config, marathon_west_files, self.cur_iter)
        self.ss_dist = ss_dist
        self.model = model

        log.debug(f'Steady-state distribution: {ss_dist}')
        log.info(f"Target steady-state flux is {ss_flux}")

        # Obtain cluster-structures
        log.debug("Obtaining cluster-structures")
        model.update_cluster_structures()

        # TODO: Do this with pathlib
        struct_directory = f"{restart_directory}/structs"
        if not os.path.exists(struct_directory):
            os.makedirs(struct_directory)

        flux_filename = f"{restart_directory}/JtargetSS.txt"
        with open(flux_filename, 'w') as fp:

            log.debug(f"Writing flux to {flux_filename}")
            fp.write(str(model.JtargetSS))
            fp.close()

        ss_filename = f"{restart_directory}/pSS.txt"
        with open(ss_filename, 'w') as fp:

            log.debug(f"Writing pSS to {ss_filename}")
            np.savetxt(fp, model.pSS)
            fp.close()

        # If this is the last run of the last restart, do nothing and exit.
        # if restart_state['runs_completed'] >= self.n_runs and restart_state['restarts_completed'] >= self.n_restarts:
        #     log.info("All restarts completed!")
        #     return

        # Construct start-state file with all structures and their weights
        # TODO: Don't explicitly write EVERY structure to disk, or this will be a nightmare for large runs.
        # However, for now, it's fine...
        log.debug("Writing structures")
        # TODO: Include start states from previous runs
        sstates_filename = f"{restart_directory}/startstates.txt"
        with open(sstates_filename, 'w') as fp:

            # Track the total number of segments iterated over
            seg_idx = 0

            log.info(f"Obtaining potential start structures ({len(model.cluster_structures.items())} avail)")

            # Loop over each set of (bin index, all the structures in that bin)
            for (msm_bin_idx, structures) in tqdm.tqdm(model.cluster_structures.items()):

                # The per-segment bin probability.
                # Map a cluster number onto a cluster INDEX, because after cleaning the cluster numbers may no longer
                # be consecutive.
                bin_prob = self.ss_dist[model.cluster_mapping[msm_bin_idx]] / len(structures)

                if bin_prob == 0:
                    log.info(f"MSM-Bin {msm_bin_idx}  has probability 0, so not saving any structs from it.")
                    continue

                # The total amount of WE weight in this MSM microbin
                msm_bin_we_weight = sum(model.cluster_structure_weights[msm_bin_idx])

                # Write each structure to disk. Loop over each structure within a bin.
                for struct_idx, structure in enumerate(structures):

                    structure_filename = (
                        f"{struct_directory}/bin{msm_bin_idx}_" f"struct{struct_idx}.{STRUCT_EXTENSIONS[self.struct_filetype]}"
                    )

                    with self.struct_filetype(structure_filename, 'w') as struct_file:

                        # One structure per segment
                        seg_we_weight = model.cluster_structure_weights[msm_bin_idx][struct_idx]

                        # Structure weights are set according to Algorithm 5.3 in
                        # Aristoff, D. & Zuckerman, D. M. Optimizing Weighted Ensemble Sampling of Steady States.
                        # Multiscale Model Sim 18, 646–673 (2020).
                        structure_weight = seg_we_weight * (bin_prob / msm_bin_we_weight)

                        topology = model.reference_structure.topology

                        try:
                            angles = model.reference_structure.unitcell_angles[0]
                            lengths = model.reference_structure.unitcell_lengths[0] * 10
                        # This throws typeerror if reference_structure.unitcell_angles is None, or AttributeError
                        #   if reference_structure.unitcell_angles doesn't exist.
                        except (TypeError, AttributeError):
                            angles, lengths = None, None

                        coords = structure * 10  # Correct units

                        # Write the structure file
                        if self.struct_filetype is md.formats.PDBTrajectoryFile:
                            struct_file.write(coords, topology, modelIndex=1, unitcell_angles=angles, unitcell_lengths=lengths)

                        elif self.struct_filetype is md.formats.AmberRestartFile:
                            # AmberRestartFile takes slightly differently named keyword args
                            struct_file.write(coords, time=None, cell_angles=angles, cell_lengths=lengths)

                        else:
                            # Otherwise, YOLO just hope all the positional arguments are in the right place
                            log.warning(
                                f"This output filetype ({self.struct_filetype}) is probably supported, "
                                f"but not explicitly handled."
                                " You should ensure that it takes argument as (coords, topology)"
                            )
                            struct_file.write(coords, topology)
                            raise Exception("Don't know what extension to use for this filetype")

                        # Add this start-state to the start-states file
                        # This path is relative to WEST_SIM_ROOT
                        fp.write(f'b{msm_bin_idx}_s{struct_idx} {structure_weight} {structure_filename}\n')
                        seg_idx += 1

        ### Start the new simulation

        bstates_str = ""
        for original_bstate in original_bstates:
            orig_bstate_prob = original_bstate.probability
            orig_bstate_label = original_bstate.label
            orig_bstate_aux = original_bstate.auxref

            bstate_str = f"{orig_bstate_label} {orig_bstate_prob} {orig_bstate_aux}\n"

            bstates_str += bstate_str

        bstates_filename = f"{restart_directory}/basisstates.txt"
        with open(bstates_filename, 'w') as fp:
            fp.write(bstates_str)

        tstates_str = ""
        for original_tstate in original_tstates:
            orig_tstate_label = original_tstate.label
            # TODO: Handle multidimensional pcoords
            orig_tstate_pcoord = original_tstate.pcoord[0]

            tstate_str = f"{orig_tstate_label} {orig_tstate_pcoord}\n"
            tstates_str += tstate_str
        tstates_filename = f"{restart_directory}/targetstates.txt"
        with open(tstates_filename, 'w') as fp:
            fp.write(tstates_str)

        # Pickle the model
        objFile = f"{restart_directory}/hamsm.obj"
        with open(objFile, "wb") as objFileHandler:
            log.debug("Pickling model")
            pickle.dump(model, objFileHandler, protocol=4)
            objFileHandler.close()

        if last_restart:
            log.info("All restarts completed! Finished.")
            return

        # Update restart_file file
        restart_state['restarts_completed'] += 1
        # If we're doing a restart, then reset the number of completed runs to 0 for the next marathon.
        restart_state['runs_completed'] = 0
        with open(self.restart_file, 'w') as fp:
            json.dump(restart_state, fp)

        log.info("Initializing new run")

        # TODO: Read this from config if available
        segs_per_state = 1

        old_initialization_path = self.initialization_file
        new_initialization_path = f"{restart_directory}/{self.initialization_file}"
        log.debug(f"Moving initialization file from {old_initialization_path} to {new_initialization_path}.")
        shutil.move(old_initialization_path, new_initialization_path)

        initialization_state = {
            'tstate_file': tstates_filename,
            'bstate_file': bstates_filename,
            'sstate_file': sstates_filename,
            'tstates': None,
            'bstates': None,
            'sstates': None,
            'segs_per_state': segs_per_state,
        }

        with open(self.initialization_file, 'w') as fp:
            json.dump(initialization_state, fp)

        westpa.rc.pstatus(
            f"\n\n"
            f"===== Restart {restart_state['restarts_completed']}, "
            + f"Run {restart_state['runs_completed']+1} initializing =====\n"
        )

        westpa.rc.pstatus(
            f"\nRun: \n\t w_init --tstate-file {tstates_filename} "
            + f"--bstate-file {bstates_filename} --sstate-file {sstates_filename} --segs-per-state {segs_per_state}\n"
        )
        w_init.initialize(**initialization_state, shotgun=False)

        log.info("New WE run ready!")
        westpa.rc.pstatus(f"\n\n===== Restart {restart_state['restarts_completed']} running =====\n")

        w_run.run_simulation()
