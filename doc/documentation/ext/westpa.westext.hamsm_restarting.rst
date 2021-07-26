westpa.westext.hamsm_restarting package
===================================

Description
-----------

This plugin leverages haMSM analysis `[1]`_ to provide simulation post-analysis. This post-analysis can be used on its own,
or can be used to initialize and run new WESTPA simulations using structures in the haMSM's best estimate of steady-state
as described in `[2]`_, which may accelerate convergence to steady-state.

haMSM analysis is performed using the `msm_we <https://github.com/jdrusso/msm_we>`_ library.

Usage
-----

Configuration
*************

:code:`west.cfg`
++++++++++++++++

This plugin requires the following section in :code:`west.cfg` (or whatever your WE configuration file is named):

.. code-block:: yaml

    plugins:
    - plugin: westpa.westext.hamsm_restarting.restart_driver.RestartDriver
      n_restarts: 0 	# Number of restarts to perform
      n_runs: 5			# Number of runs within each restart
      coord_len: 2		# Length of pcoords returned
      initialization_file: restart_initialization.json	# JSON describing w_run parameters for new runs
      ref_pdb_file: common_files/bstate.pdb 		    # File containing reference structure
      init_pdb_file: bstates/bstate.pdb 	            # File containing basis state reference
      model_name: NaClFlux		    # Name for msm_we model
      n_clusters: 25				# Number of clusters in haMSM building
      we_folder: .					# Should point to the same directory as WEST_SIM_ROOT
      target_pcoord1: 2.6			# Progress value in the target state
      basis_pcoord1_min: 14.0		# Maximum progress coordinate within the basis state
      basis_pcoord1_max: 15.0		# Minimum progress coordinate within the basis state
      pcoord_ndim0: 1	            # Dimensionality of progress coordinate
      dim_reduce_method: pca		# Dimensionality reduction scheme, either "pca", "vamp", or "none"
      parent_traj_filename: parent.xml		# Name of parent file in each segment
      child_traj_filename: seg.xml		    # Name of child file in each segment
      user_functions: westpa_scripts/restart_overrides.py	# Python file defining coordinate processing
      struct_filetype: mdtraj.formats.PDBTrajectoryFile 	# Filetype for output start-structures

Some sample parameters are provided in the above, but of course should be modified to your specific system.

Note that ref_pdb_file and init_pdb_file can be any filetype supported by :code:`msm_we.initialize()`'s structure loading.
At the time of writing, this is limited to PDB, however that is planned to be extended.

Also in this file, :code:`west.data.data_refs.basis_state` MUST point to
:code:`$WEST_SIM_ROOT/{basis_state.auxref}` and not a subdirectory if restarts are being used.
This is because when the plugin initiates a restart, start_state
references in :code:`$WEST_SIM_ROOT/restartXX/start_states.txt` are set relative to :code:`$WEST_SIM_ROOT`. All basis/start
state references are defined relative to :code:`west.data.data_refs.basis_state`, so if that points to a subdirectory of
:code:`$WEST_SIM_ROOT`, those paths will not be accurate.

:code:`restart_initialization.json`
+++++++++++++++++++++++++++++++++++

.. code-block:: json

    {
        "bstates":["start,1,bstates/bstate.pdb"],
        "tstates":["bound,2.6"],
        "bstate-file":"bstates/bstates.txt",
        "tstate-file" :"tstate.file",
        "segs-per-state": 1
    }

It is not necessary to specify both in-line states and a state-file for each, but that is shown in the sample for
completeness.

It is important that :code:`bstates` and :code:`tstates` are lists of strings, and not just strings, even if only one
bstate/tstate is being used!

With :code:`n_runs > 1`, before doing any restart, multiple independent runs are performed. However, before the first
restart (this applies if no restarts are performed as well), the plugin has no way of accessing the parameters that
were initially passed to :code:`w_init` and :code:`w_run`.

Therefore, it is necessary to store those parameters in a file, so the plugin can read them and initiate subsequent runs.

After the first restart is performed, the plugin writes this file itself, so it is only necessary to manually configure
for that first set of runs.

Featurization overrides
+++++++++++++++++++++++

.. code-block:: python


    import numpy as np
    import mdtraj as md

    def processCoordinates(self, coords):
            log.debug("Processing coordinates")

            if self.dimReduceMethod == "none":
                nC = np.shape(coords)
                nC = nC[0]
                ndim = 3 * self.nAtoms
                data = coords.reshape(nC, 3 * self.nAtoms)
                return data

            if self.dimReduceMethod == "pca" or self.dimReduceMethod == "vamp":

                ### NaCl RMSD dimensionality reduction
                log.warning("Hardcoded selection: Doing dim reduction for Na, Cl. This is only for testing!")
                indNA = self.reference_structure.topology.select("element Na")
                indCL = self.reference_structure.topology.select("element Cl")

                diff = np.subtract(coords[:, indNA], coords[:, indCL])

                dist = np.array(np.sqrt(
                    np.mean(
                        np.power(
                            diff,
                            2)
                    , axis=-1)
                ))

                return dist

This is the file whose path is provided in the configuration file in :code:`plugin.user_functions`, and must be a Python
file defining a function named :code:`processCoordinates(self, coords)` which takes a numpy array of coordinates,
featurizes it, and returns the numpy array of feature-coordinates.

This is left to be user-provided because whatever featurization you do will be system-specific. The provided function
is monkey-patched into the :code:`msm_we.modelWE` class.

An example is provided above, which does a simple RMSD coordinate reduction for the NaCl association tutorial system.

Doing only post-analysis
************************

If you want to ONLY use this for haMSM post-analysis, and not restarting, just set :code:`n_restarts: 0` in the configuration.


Work manager for restarting
***************************

If you're using some parallelism (which you should), and you're using the plugin to do restarts or multiple runs,
then your choice of work manager can be important.
This plugin handles starting new WESTPA runs using the Python API.
The process work manager, by default, uses :code:`fork` to start new workers which seems to eventually causes
memory issues, since :code:`fork` passes the entire contents of the parent to each child.
Switching the spawn method to :code:`forkserver` or :code:`spawn` may introduce other issues.

Using the ZMQ work manager works well. The MPI work manager should also work well, though is untested.
Both of these handle starting new workers in a more efficient way, without copying the full state of the parent.

Potential Pitfalls/Troubleshooting
******************
- Verify that `msm_we <https://github.com/jdrusso/msm_we>`_ is installed

- Verify that :code:`restart_initialization.json` has been correctly set

- This plugin does not yet attempt to resolve environment variables in the config, so things like say, $WEST_SIM_ROOT, will be interpreted literally in paths

References
**********
_`[1]` Suárez, E., Adelman, J. L. & Zuckerman, D. M. Accurate Estimation of Protein Folding and Unfolding Times: Beyond Markov State Models. J Chem Theory Comput 12, 3473–3481 (2016).

_`[2]` Copperman, J. & Zuckerman, D. M. Accelerated Estimation of Long-Timescale Kinetics from Weighted Ensemble Simulation via Non-Markovian “Microbin” Analysis. J Chem Theory Comput 16, 6763–6775 (2020).


API
---

westpa.westext.hamsm\_restarting module
***************************************

.. automodule:: westpa.westext.hamsm_restarting.restart_driver
   :members:
   :undoc-members:
   :show-inheritance:
