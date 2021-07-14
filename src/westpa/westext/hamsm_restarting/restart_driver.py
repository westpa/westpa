import logging
import operator
import numpy as np

import westpa
from westpa.cli.core import w_init

import json

import os
import subprocess
import shutil
import pickle

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

FORMAT = "%(message)s"
logging.basicConfig(level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
log = logging.getLogger("restart_driver")
log.setLevel("DEBUG")
logging.getLogger("msm_we").setLevel("INFO")


# TODO: Break this out into a separate module, let it be specified (if it's necessary) as a plugin option
#   This may not always be required -- i.e. you may be able to directly output to the h5 file in your propagator
def prepare_coordinates(plugin_config, h5file):
    """
    Copy relevant coordinates from trajectory files into <iteration>/auxdata/coord of the h5 file.

    Directly modifies the input h5 file.

    Adapted from original msmWE collectCoordinates.py script.

    Parameters
    ----------
    plugin_config: YAMLConfig object
        Stores the configuration options provided to the plugin in the WESTPA configuration file

    h5file: h5py.File
        WESTPA h5 data file
    """

    fileSpecifier = plugin_config.get('file_specifier')
    refPDBfile = plugin_config.get('ref_pdb_file')
    initPDBfile = plugin_config.get('init_pdb_file')
    modelName = plugin_config.get('model_name')

    # TODO: Don't need this explicit option, use WEST_SIM_ROOT or something
    WEfolder = plugin_config.get('we_folder')

    parentTraj = plugin_config.get('parent_traj_filename')
    childTraj = plugin_config.get('child_traj_filename')

    model = msm_we.modelWE()
    log.info('Preparing coordinates...')
    model.initialize(fileSpecifier, refPDBfile, initPDBfile, modelName)
    model.get_iterations()

    for n_iter in tqdm.tqdm(range(1, model.maxIter)):

        nS = model.numSegments[n_iter - 1].astype(int)
        coords = np.zeros((nS, 2, model.nAtoms, 3))
        dsetName = "/iterations/iter_%08d/auxdata/coord" % int(n_iter)

        coords_exist = False
        try:
            dset = h5file.create_dataset(dsetName, np.shape(coords))
        except RuntimeError:
            log.debug('coords exist for iteration ' + str(n_iter) + ' NOT overwritten')
            coords_exist = True
            continue

        for iS in range(nS):
            trajpath = WEfolder + "/traj_segs/%06d/%06d" % (n_iter, iS)

            # TODO: HACK, this should  not be hardcoded
            try:
                coord0 = np.squeeze(md.load(f'{trajpath}/{parentTraj}', top=model.reference_structure.topology)._xyz)
            except OSError:
                # coord0 = np.squeeze(md.load(trajpath + '/parent.xml', top=model.reference_structure.topology)._xyz)
                log.warning("Parent traj file doesn't exist, loading reference structure coords")
                coord0 = np.squeeze(model.reference_structure._xyz)

            coord1 = np.squeeze(md.load(f'{trajpath}/{childTraj}', top=model.reference_structure.topology)._xyz)

            coords[iS, 0, :, :] = coord0
            coords[iS, 1, :, :] = coord1

        if not coords_exist:
            # log.debug(f"Writing coords for iter {n_iter}")
            dset[:] = coords


def msmwe_compute_ss(plugin_config, last_iter):
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
    model = msm_we.modelWE()

    fileSpecifier = plugin_config.get('file_specifier')
    refPDBfile = plugin_config.get('ref_pdb_file')
    initPDBfile = plugin_config.get('init_pdb_file')
    modelName = plugin_config.get('model_name')
    n_clusters = plugin_config.get('n_clusters')
    # last_iter = plugin_config.get('last_iter')

    # Fire up the model object
    # (Eventually this will just go in __init__)
    model.initialize(fileSpecifier, refPDBfile, initPDBfile, modelName)

    # Set some model parameters
    model.WEtargetp1 = plugin_config.get('target_pcoord1')
    model.WEbasisp1_min = plugin_config.get('basis_pcoord1_min')
    model.WEbasisp1_max = plugin_config.get('basis_pcoord1_max')
    model.pcoord_ndim0 = plugin_config.get('pcoord_ndim0')
    model.dimReduceMethod = plugin_config.get('dim_reduce_method')

    model.n_lag = n_lag

    last_iter_cluster = last_iter  # model.maxIter-1 #last iter often not complete
    #####
    # TODO: magic number? Maximum possible number of segments I think based on the while loop below
    i = last_iter_cluster
    coordSet = np.zeros((0, model.nAtoms, 3))  # extract coordinate libraries for clustering

    log.debug("Loading in iteration data.. (this could take a while)")
    log.debug(f'coord shape is {coordSet.shape}')

    ##### Replacement
    # First dimension is the total number of segments
    model.get_iterations()
    # Ignore the last iteration, as above
    total_segments = int(sum(model.numSegments[:-1]))
    coordSet = np.zeros((total_segments, model.nAtoms, 3))  # extract coordinate libraries for clustering
    pcoordSet = np.zeros((total_segments, model.pcoord_ndim))

    last_seg = total_segments

    # Update iterations N+1 -> 1
    for i in tqdm.tqdm(range(last_iter, 0, -1)):
        # log.debug(f"Appending coords from iteration {i}/{last_iter} (goes backwards to 0)")

        model.load_iter_data(i)
        model.get_iter_coordinates()

        first_seg = last_seg - len(model.segindList)
        assert first_seg >= 0, "Referencing a segment that doesn't exist"

        # log.debug(f"This covers segments {first_seg} to {last_seg}")

        indGood = np.squeeze(np.where(np.sum(np.sum(model.cur_iter_coords, 2), 1) != 0))

        coordSet[first_seg:last_seg] = model.cur_iter_coords[indGood, :, :]
        pcoordSet[first_seg:last_seg] = model.pcoord1List[indGood, :]

        last_seg = first_seg

    log.debug(f"Segment weights has {list(model.seg_weights.keys())}")
    log.debug(f"Segment weights has {len(list(model.seg_weights.keys()))}")
    # raise Exception
    #####

    # Set the coords, and pcoords
    # TODO: There's no need to set these as attributes of the object
    model.all_coords = coordSet
    model.pcoordSet = pcoordSet

    log.debug(f"CoordLIST: {model.cur_iter_coords.shape}")
    log.debug(f"CoordSET: {model.all_coords.shape}")

    # TODO: Are first_iter and last_iter used consistently everywhere? Some places they're taken as parameters,
    #   some places the current value is just pulled from state
    # TODO: What does "cluster" mean?
    first_iter_cluster = i
    model.first_iter = first_iter_cluster
    model.last_iter = last_iter_cluster

    # TODO: Related to above comment, just use coordSet not model.coordSet
    n_coords = np.shape(model.all_coords)[0]

    model.dimReduce()

    clusterFile = (
        modelName + "_clusters_s" + str(first_iter_cluster) + "_e" + str(last_iter_cluster) + "_nC" + str(n_clusters) + ".h5"
    )
    exists = os.path.isfile(clusterFile)
    exists = False
    log.warning("Skipping any potential cluster reloading!")

    # If a cluster file with the name corresponding to these parameters exists, load clusters from it.
    if exists:
        print("loading clusters...")
        model.load_clusters(clusterFile)
    # Otherwise, do the clustering (which will create and save to that file)
    else:
        # FIXME: This gives the wrong shape, but loading from the clusterfile gives the right shape
        print("clustering " + str(n_coords) + " coordinates into " + str(n_clusters) + " clusters...")
        model.cluster_coordinates(n_clusters)

    first_iter = 1
    model.get_fluxMatrix(n_lag, first_iter, last_iter)  # extracts flux matrix, output model.fluxMatrixRaw
    westpa.rc.pstatus(f"Unprocessed flux matrix has shape {model.fluxMatrixRaw.shape}")
    model.organize_fluxMatrix()  # gets rid of bins with no connectivity, sorts along p1, output model.fluxMatrix
    model.get_Tmatrix()  # normalizes fluxMatrix to transition matrix, output model.Tmatrix

    westpa.rc.pstatus("Flux matrix:")
    westpa.rc.pstatus(model.fluxMatrix)
    westpa.rc.pstatus("T matrix: ")
    westpa.rc.pstatus(model.Tmatrix)
    westpa.rc.pstatus(f"Processed flux matrix has shape {model.fluxMatrix.shape}")

    logging.getLogger("msm_we").setLevel("DEBUG")
    model.get_steady_state_algebraic()  # gets steady-state from eigen decomp, output model.pSS
    logging.getLogger("msm_we").setLevel("INFO")
    model.get_steady_state_target_flux()  # gets steady-state target flux, output model.JtargetSS

    # FIXME: Why is this the wrong shape?
    ss_alg = model.pSS
    ss_flux = model.JtargetSS

    westpa.rc.pstatus("Got steady state:")
    westpa.rc.pstatus(ss_alg)
    westpa.rc.pstatus(ss_flux)

    westpa.rc.pstatus("Completed flux matrix calculation and steady-state estimation!")

    return ss_alg, ss_flux, model


def reduce_array(Aij):
    """
    Remove empty rows and columns from an array Aij and return the reduced
        array Bij and the list of non-empty states,

    Parameters
    ----------
    Aij: np.ndarray
        The array to reduce.

    Returns
    -------
    Bij: np.ndarray
        Aij, with empty rows and columns removed

    nonempty: list
        Indices of the empty rows and columns which were removed from Aij
    """

    nonempty = list(range(0, Aij.shape[0]))
    eps = np.finfo(Aij.dtype).eps

    for i in range(0, Aij.shape[0]):
        if (Aij[i, :] < eps).all() and (Aij[:, i] < eps).all():
            nonempty.pop(nonempty.index(i))

    nne = len(nonempty)
    Bij = np.zeros((nne, nne))

    for i in range(0, nne):
        for j in range(0, nne):
            Bij[i, j] = Aij[nonempty[i], nonempty[j]]

    return Bij, nonempty


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

        self.coord_len = plugin_config.get('coord_len', 2)
        self.n_restarts = plugin_config.get('n_restarts', 2)

        # This should be low priority, because it closes the H5 file and starts a new WE run. So it should run LAST
        #   after any other plugins.
        self.priority = plugin_config.get('priority', 100)  # I think a big number is lower priority...

        # sim_manager.register_callback(sim_manager.post_we,
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

        The algorithm is as follows:
            1. Check to see if we've just completed the final iteration
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

        restart_state = {'restarts_completed': 0}

        # Look for a restart.dat file to get the current state (how many restarts have been performed already)
        if os.path.exists(self.restart_file):
            with open(self.restart_file, 'r') as fp:
                restart_state = json.load(fp)

            # If we just completed the simulation of the final restart, do nothing and exit, we're all done
            if restart_state['restarts_completed'] >= self.n_restarts:
                log.info("All restarts completed! Exiting.")
                return

        log.debug(f"{restart_state['restarts_completed']}/{self.n_restarts} completed")

        # Build the haMSM
        westpa.rc.pstatus("Initializing haMSM")
        prepare_coordinates(self.plugin_config, self.data_manager.we_h5file)

        # Get haMSM steady-state estimate
        westpa.rc.pstatus("Computing steady-state")
        ss_dist, ss_flux, model = msmwe_compute_ss(self.plugin_config, self.cur_iter)
        self.ss_dist = ss_dist

        # For some reason, this is returned irregularly as a list of lists sometimes.
        # TODO: Fix that...
        westpa.rc.pstatus(f"Flattening {ss_dist}")
        # self.ss_dist = np.array(list(flatten(ss_dist)))
        self.model = model

        # Obtain cluster-structures
        westpa.rc.pstatus("Obtaining cluster-structures")

        # TODO: I think when the flux matrix is cleaned, the cluster structures are not
        #  reassigned to the new, reduced set of clusters
        # logging.getLogger("msm_we").setLevel("DEBUG")
        model.update_cluster_structures()
        # logging.getLogger("msm_we").setLevel("INFO")

        # Construct start-state file with all structures and their weights
        # TODO: Don't explicitly write EVERY structure to disk, or this will be a nightmare for large runs.
        # However, for now, it's fine...
        log.debug("Writing structures")

        # TODO: Do this with pathlib
        restart_directory = f"restart{restart_state['restarts_completed']}"
        struct_directory = f"{restart_directory}/structs"
        if not os.path.exists(struct_directory):
            os.makedirs(struct_directory)

        log.debug(f'Steady-state distribution: {self.ss_dist}')

        log.info("Computing fluxes")
        model.get_steady_state_algebraic()  # gets steady-state from eigen decomp, output model.pSS
        model.get_steady_state_target_flux()  # gets steady-state target flux, output model.JtargetSS

        log.info(f"Target steady-state flux is {model.JtargetSS}")

        flux_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_JtargetSS.txt"
        with open(flux_filename, 'w') as fp:

            log.info(f"Writing flux to {flux_filename}")
            # np.savetxt(fp, model.JtargetSS)
            fp.write(str(model.JtargetSS))
            fp.close()

        ss_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_pSS.txt"
        with open(ss_filename, 'w') as fp:

            log.info(f"Writing pSS to {ss_filename}")
            # fp.write(model.pSS)
            np.savetxt(fp, model.pSS)
            fp.close()

        # TODO: Include start states from previous runs
        sstates_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_startstates.txt"
        with open(sstates_filename, 'w') as fp:

            # Track the total number of segments iterated over
            seg_idx = 0

            # This comparison is not correct
            # n_structures = sum([len(x) for x in model.cluster_structures.values()])
            # # Make sure the total number of structures matches the number of structure weights
            # #   Do this here, because if we do it above we have to get the list of structures and then flatten it.
            # log.debug(f"{n_structures} structures and {len(model.seg_weights)} weights")
            # assert n_structures == len(model.cluster_structure_weights), "Number of structures doesn't match number of weight entries"

            log.info("Obtaining potential start structures")

            # Loop over each set of (bin index, all the structures in that bin)
            for (msm_bin_idx, structures) in tqdm.tqdm(model.cluster_structures.items()):

                # The per-segment bin probability
                bin_prob = self.ss_dist[msm_bin_idx] / len(structures)

                if bin_prob == 0:
                    westpa.rc.pstatus(f"MSM-Bin {msm_bin_idx}  has probability 0, so not saving any structs from it.")
                    continue

                msm_bin_we_weight = sum(model.cluster_structure_weights[msm_bin_idx])

                # Write each structure to disk
                # TODO: I need atomtypes to go with each coordinate!
                #   Do I? I don't think that's true.
                # Loop over each structure  within a bin.
                for struct_idx, structure in enumerate(structures):

                    structure_filename = f"{struct_directory}/bin{msm_bin_idx}_struct{struct_idx}.pdb"
                    with md.formats.PDBTrajectoryFile(structure_filename, 'w') as struct_file:

                        # One structure per segment
                        seg_we_weight = model.cluster_structure_weights[msm_bin_idx][struct_idx]

                        structure_weight = bin_prob * (seg_we_weight / msm_bin_we_weight)

                        topology = model.reference_structure.topology
                        angles = model.reference_structure.unitcell_angles[0]
                        lengths = model.reference_structure.unitcell_lengths[0] * 10
                        # log.debug(angles)
                        # log.debug(lengths)
                        coords = structure * 10  # Correct units

                        # Write the structure file
                        struct_file.write(coords, topology, modelIndex=1, unitcell_angles=angles, unitcell_lengths=lengths)

                        # Add this start-state to the start-states file
                        # This path is relative to WEST_SIM_ROOT
                        fp.write(f'b{msm_bin_idx}_s{struct_idx} {structure_weight} {structure_filename}\n')
                        seg_idx += 1

        ### Start the new simulation
        # Back up west.h5
        shutil.copyfile('west.h5', f"{restart_directory}/west_restart{restart_state['restarts_completed']}.h5")

        # Keep west.cfg as-is, it should be reusable in the next restart w/o modification!

        ## Run init w/ the new start-state file + the original basis state
        # I can get these from sim_manager.current_iter_bstates
        original_bstates = self.sim_manager.current_iter_bstates

        if original_bstates is None:
            original_bstates = self.data_manager.get_basis_states(self.sim_manager.n_iter - 1)

        assert original_bstates is not None, "Bstates are none in the current iteration"

        bstates_str = ""
        for original_bstate in original_bstates:
            orig_bstate_prob = original_bstate.probability
            orig_bstate_label = original_bstate.label
            orig_bstate_aux = original_bstate.auxref

            bstate_str = f"{orig_bstate_label} {orig_bstate_prob} {orig_bstate_aux}\n"

            bstates_str += bstate_str

        bstates_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_basisstates.txt"
        with open(bstates_filename, 'w') as fp:
            fp.write(bstates_str)

        # TODO: Don't repeat this code!
        original_tstates = self.data_manager.get_target_states(self.cur_iter)
        tstates_str = ""
        for original_tstate in original_tstates:
            orig_tstate_label = original_tstate.label
            # TODO: Handle multidimensional pcoords
            orig_tstate_pcoord = original_tstate.pcoord[0]

            tstate_str = f"{orig_tstate_label} {orig_tstate_pcoord}\n"
            tstates_str += tstate_str
        tstates_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_targetstates.txt"
        with open(tstates_filename, 'w') as fp:
            fp.write(tstates_str)

        # Close west.h5
        # TODO: Doing this here is kinda janky. But it shouldn't really be a problem, especially as long as
        #   this plugin's priority is low...
        # data_manager.finalize_run() currently only flushes and closes the west.h5, and at time of writing,
        # nothing else happens after that, WESTPA just finishes running.
        # But you could imagine a situation where something *is* supposed to write to west.h5 in finalize step.
        # In that case, if this has higher priority, the other thing would crash.
        self.data_manager.finalize_run()

        westpa.rc.pstatus(
            f"Run: \n\t w_init --tstate-file {tstates_filename} "
            + f"--bstate-file {bstates_filename} --sstate-file {sstates_filename} --segs-per-state 5"
        )

        # TODO: Move traj_segs and seg_logs

        for data_folder in ['traj_segs', 'seg_logs']:
            old_path = data_folder
            new_path = f"{restart_directory}/{data_folder}"

            log.debug(f"Moving {old_path} to {new_path}")

            try:
                shutil.move(old_path, new_path)
            except FileNotFoundError:
                log.warning(f"Folder {data_folder}  was not found. This may be normal, but check your configuration.")
                continue
            else:
                os.mkdir(old_path)

        log.info("Initializing new run")

        westpa.rc.pstatus(f"\n\n===== Restart {restart_state['restarts_completed']} initializing =====\n")
        w_init.initialize(
            tstate_file=tstates_filename,
            bstate_file=bstates_filename,
            sstate_file=sstates_filename,
            tstates=None,
            bstates=None,
            sstates=None,
            segs_per_state=5,
            shotgun=False,
        )

        # Update restart_file file
        restart_state['restarts_completed'] += 1
        with open(self.restart_file, 'w') as fp:
            json.dump(restart_state, fp)

        # Pickle the model
        objFile = f"{restart_directory}/restart{restart_state['restarts_completed']}_hamsm.obj"
        with open(objFile, "wb") as objFileHandler:
            del model.clusters
            log.debug("Pickling model")
            pickle.dump(model, objFileHandler)
            objFileHandler.close()

        log.info("New WE run ready!")
        westpa.rc.pstatus(f"\n\n===== Restart {restart_state['restarts_completed']} running =====\n")

        # TODO: Do this via the Python API instead of by running a run.sh
        #       Get the flags that were passed to w_run for this one, and pass it to the next.
        current_env = os.environ.copy()
        subprocess.Popen('./run.sh', env=current_env).wait()
