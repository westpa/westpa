import logging
import operator
import numpy as np

import westpa

import json

import os
import shutil
import sys

sys.path.append("/Users/russojd/Research/molssi_2021/msmWE_cc/msm_we")
from msm_we import msm_we

import mdtraj as md

EPS = np.finfo(np.float64).eps

from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
log = logging.getLogger("restart_driver")
log.setLevel("DEBUG")
logging.getLogger("msm_we").setLevel("INFO")

from collections.abc import Iterable


def flatten(_list):
    for el in _list:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


# This is basically CollectCoordinates
# This is an extra step that, if necessary, puts all
# TODO: Break this out into a separate module, let it be specified (if it's necessary) as a plugin option
#   This may not always be required -- i.e. you may be able to directly output to the h5 file in your propagator
def prepare_coordinates(plugin_config, h5file):
    """
    Copy relevant coordinates from trajectory files into auxdata/coord of the h5 file
    """

    fileSpecifier = plugin_config.get('file_specifier')
    refPDBfile = plugin_config.get('ref_pdb_file')
    initPDBfile = plugin_config.get('init_pdb_file')
    modelName = plugin_config.get('model_name')

    # TODO: Don't need this explicit option, use WEST_SIM_ROOT or something
    WEfolder = plugin_config.get('we_folder')

    model = msm_we.modelWE()
    log.debug('initializing...')
    model.initialize(fileSpecifier, refPDBfile, initPDBfile, modelName)
    model.get_iterations()

    for n_iter in range(1, model.maxIter):

        nS = model.numSegments[n_iter - 1].astype(int)
        coords = np.zeros((nS, 2, model.nAtoms, 3))
        for iS in range(nS):
            trajpath = WEfolder + "/traj_segs/%06d/%06d" % (n_iter, iS)
            coord0 = np.squeeze(md.load(trajpath + '/parent.xml', top=model.reference_structure.topology)._xyz)
            coord1 = np.squeeze(md.load(trajpath + '/seg.xml', top=model.reference_structure.topology)._xyz)
            coords[iS, 0, :, :] = coord0
            coords[iS, 1, :, :] = coord1

        dsetName = "/iterations/iter_%08d/auxdata/coord" % int(n_iter)
        try:
            dset = h5file.create_dataset(dsetName, np.shape(coords))
            dset[:] = coords
        except Exception:
            log.debug('coords exist for iteration ' + str(n_iter) + ' NOT overwritten\n')

    # f.close()


def msmwe_compute_ss(plugin_config, last_iter):

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
    # TODO: magic number? Maximum possible number of segments I think based on the while loop below
    ncoords = 100000
    i = last_iter_cluster

    numCoords = 0
    coordSet = np.zeros((0, model.nAtoms, 3))  # extract coordinate libraries for clustering
    pcoordSet = np.zeros((0, model.pcoord_ndim))

    log.debug("Loading in iteration data.. (this could take a while)")

    # Keep adding iterations until the total number of segments is equal to ncoords
    while numCoords < ncoords:
        print(f"Doing iteration {i}/{last_iter_cluster} (goes backwards to 0) | {numCoords}", end='\r')
        # TODO: Replace this, which iterates over a set of trajectory files, with something that just looks at auxdata
        # Load data into the model's state
        model.load_iter_data(i)
        model.get_iter_coordinates()

        # TODO: Slowwwwwww appends
        indGood = np.squeeze(np.where(np.sum(np.sum(model.cur_iter_coords, 2), 1) != 0))
        coordSet = np.append(coordSet, model.cur_iter_coords[indGood, :, :], axis=0)
        pcoordSet = np.append(pcoordSet, model.pcoord1List[indGood, :], axis=0)
        numCoords = np.shape(coordSet)[0]
        i = i - 1

        if i == 0:
            break

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

    model.get_steady_state_algebraic()  # gets steady-state from eigen decomp, output model.pSS
    model.get_steady_state_target_flux()  # gets steady-state target flux, output model.JtargetSS

    ss_alg = model.pSS
    ss_flux = model.JtargetSS

    westpa.rc.pstatus("Got steady state:")
    westpa.rc.pstatus(ss_alg)
    westpa.rc.pstatus(ss_flux)

    # objFile = (
    #         modelName
    #         + "_s"
    #         + str(first_iter)
    #         + "_e"
    #         + str(last_iter)
    #         + "_nC"
    #         + str(n_clusters)
    #         + ".obj"
    # )
    # objFileHandler = open(objFile, "wb")
    # del model.clusters
    # # print("Pickling model")
    # pickle.dump(model, objFileHandler)
    # objFileHandler.close()

    westpa.rc.pstatus("Completed flux matrix calculation and steady-state estimation!")

    return ss_alg, ss_flux, model


def reduce_array(Aij):
    """Remove empty rows and columns from an array Aij and return the reduced
        array Bij and the list of non-empty states"""

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
    def __init__(self, sim_manager, plugin_config):

        westpa.rc.pstatus("Restart plugin initialized")

        if not sim_manager.work_manager.is_master:
            westpa.rc.pstatus("Reweighting not master, skipping")
            return

        self.data_manager = sim_manager.data_manager
        self.sim_manager = sim_manager

        self.plugin_config = plugin_config

        self.restart_file = plugin_config.get('restart_file', 'restart.dat')

        self.coord_len = plugin_config.get('coord_len', 2)
        # self.reweight_period = plugin_config.get('reweight_period', 2)
        self.n_restarts = plugin_config.get('n_restarts', 2)
        self.priority = plugin_config.get('priority', 100)  # I think a big number is lower priority...

        # sim_manager.register_callback(sim_manager.post_we,
        sim_manager.register_callback(sim_manager.finalize_run, self.prepare_new_we, self.priority)

        # Initialize data
        self.ss_alg = None
        self.ss_dist = None
        self.model = None

    def get_original_bins(self):
        """
        Returns the WE bins and their probabilities at the end of the previous iteration.
        """

        we_driver = self.sim_manager.we_driver
        bins = we_driver.next_iter_binning
        n_bins = len(bins)
        binprobs = np.fromiter(map(operator.attrgetter('weight'), bins), dtype=np.float64, count=n_bins)

        return bins, binprobs

    @property
    def cur_iter(self):
        """
        Returns
        -------
        int: The current iteration. Subtract one, because in finalize_run the iter has been incremented
        """
        return self.sim_manager.n_iter - 1

    @property
    def is_last_iteration(self):
        """
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
                raise Exception
                return

        log.debug(f"{restart_state['restarts_completed']}/{self.n_restarts} completed")

        # Build the haMSM
        westpa.rc.pstatus("Initializing haMSM")
        prepare_coordinates(self.plugin_config, self.data_manager.we_h5file)

        # Get haMSM steady-state estimate
        westpa.rc.pstatus("Computing steady-state")
        ss_dist, ss_flux, model = msmwe_compute_ss(self.plugin_config, self.cur_iter)

        # For some reason, this is returned irregularly as a list of lists sometimes.
        # TODO: Fix that...
        westpa.rc.pstatus(f"Flattening {ss_dist}")
        self.ss_dist = np.array(list(flatten(ss_dist)))
        self.model = model

        # Obtain cluster-structures
        westpa.rc.pstatus("Obtaining cluster-structures")

        # TODO: I think when the flux matrix is cleaned, the cluster  structures are not reassigned to the new, reduced set of clusters
        logging.getLogger("msm_we").setLevel("DEBUG")
        model.update_cluster_structures()
        logging.getLogger("msm_we").setLevel("INFO")

        # Construct start-state file with all structures and their weights
        # TODO: Don't explicitly write EVERY structure to disk, or this will be a nightmare for large runs.
        # However, for now, it's fine...
        westpa.rc.pstatus("Writing structures")

        # TODO: Do this with pathlib
        restart_directory = f"restart{restart_state['restarts_completed']}"
        struct_directory = f"{restart_directory}/structs"
        if not os.path.exists(struct_directory):
            os.makedirs(struct_directory)

        log.debug('===')
        log.debug(self.ss_dist)
        log.debug('===')

        # TODO: Include start states from previous runs
        sstates_filename = f"{restart_directory}/restart{restart_state['restarts_completed']}_startstates.txt"
        with open(sstates_filename, 'w') as fp:

            for (bin_idx, structures) in model.cluster_structures.items():

                print(self.ss_dist)
                bin_prob = self.ss_dist[bin_idx] / len(structures)

                if bin_prob == 0:
                    westpa.rc.pstatus(f"MSM-Bin {bin_idx}  has probability 0, so not saving any structs from it.")
                    continue

                # Write each structure to disk
                # TODO: I need atomtypes to go with each coordinate!
                for struct_idx, structure in enumerate(structures):

                    structure_filename = f"{struct_directory}/bin{bin_idx}_struct{struct_idx}.pdb"
                    with md.formats.PDBTrajectoryFile(structure_filename, 'w') as struct_file:

                        topology = model.reference_structure.topology
                        angles = model.reference_structure.unitcell_angles[0]
                        lengths = model.reference_structure.unitcell_lengths[0] * 10
                        # log.debug(angles)
                        # log.debug(lengths)
                        coords = structure * 10  # Correct units
                        # print(coords)
                        # raise Exception

                        # Write the structure file
                        struct_file.write(coords, topology, modelIndex=1, unitcell_angles=angles, unitcell_lengths=lengths)

                        # Add this start-state to the start-states file
                        # fp.write(f'b{bin_idx}_s{struct_idx} {bin_prob} {os.path.abspath(structure_filename)}\n')
                        # This path is relative to WEST_SIM_ROOT
                        fp.write(f'b{bin_idx}_s{struct_idx} {bin_prob} {structure_filename}\n')

        ### Start the new simulation
        # Back up west.h5
        shutil.copyfile('west.h5', f"{restart_directory}/west_restart{restart_state['restarts_completed']}.h5")

        # Keep west.cfg as-is, it should be reusable in the next restart w/o modification!

        ## Run init w/ the new start-state file + the original basis state
        # I can get these from sim_manager.current_iter_bstates
        original_bstates = self.sim_manager.current_iter_bstates
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

        # with open(sstates_filename, 'a') as fp:
        #     fp.write(bstates_str)

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
        # TODO: Doing this here is kinda janky. It shouldn't really be a problem, especially as long as
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
        # w_init.initialize(tstate_file=tstates_filename,
        #                   bstate_file=bstates_filename,
        #                   sstate_file=sstates_filename,
        #                   tstates=None,
        #                   bstates=None,
        #                   sstates=None,
        #                   segs_per_state=5,
        #                   shotgun=False)

        # Update restart_file file
        restart_state['restarts_completed'] += 1
        with open(self.restart_file, 'w') as fp:
            json.dump(restart_state, fp)

        # raise Exception
