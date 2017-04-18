Checklist 
==========

Configuring a WESTPA Simulation
--------------------------------

- Files for dynamics propagation

  + Have you set up all of the files for propagating the dynamics (e.g. for
    GROMACS, the .top, .gro, .mdp, and .ndx files)?

- System implementation (``system.py``)

  + Is ``self.pcoord_len`` set to the number of data points that
    corresponds to the frequency with which the dynamics engine outputs the
    progress coordinate? Note: Many MD engines (e.g. GROMACS) output the
    initial point (i.e. zero).
  + Are the bins in the expected positions? You can easily view the positions
    of the bins using a Python interpreter.

- Initializing the simulation (``init.sh``)

  + Is the directory structure for the trajectory output files
    consistent with specifications in the master configuration file
    (``west.cfg``)?
  + Are the basis (bstate) states, and if applicable, target states (tstate),
    specified correctly?

- Calculating the progress coordinate for initial states (``get_pcoord.sh``)

  + Ensure that the procedure to extract the progress coordinate works by
    manually checking the procedure on one (or more) basis state files.
  + If your initialization (``init.sh``) gives an error message indicating the
    "incorrect shape" of the progress coordinate, check that get_pcoord.sh is
    not writing to a single file. If this is the case, w_init will crash since
    multiple threads will be simultaneously writing to a single file. To fix
    this issue, you can add $$ to the file name (e.g. change ``OUT=dist.xvg``
    to ``OUT=dist_$$.xvg``) in ``get_pcoord.sh``.

- Segment implementation (``runseg.sh``)

  + Ensure that the progress coordinate is being calculated correctly.
    If necessary, manually run a single dynamics segment (τ) for a single
    trajectory walker to do so (e.g. for GROMACS, run the .tpr file for a
    length of τ). Double check that if any analysis programs are being run
    that their input is correct.
  + Are you feeding the velocities and state information required for the
    thermostat and barostat from one dynamics segment to the next? In GROMACS,
    this information is stored in the .edr and .trr files.
    
- Log of simulation progress (``west.h5``)

  + Check that the first iteration has been initialized, i.e. typing::

      h5ls west.h5/iterations

    at the command line gives::

      iter_00000001            Group

  + In addition, the progress coordinate should be initialized as well, i.e.
    using the command::

      h5ls -d west.h5/iterations/iter_00000001/pcoord

    shows that the array is populated by zeros and the first point is the value
    calculated by get_pcoord.sh::

      pcoord                   Dataset {10, 21, 1}
         Data:
             (0,0,0) 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             (2,15,0) 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0,
             (5,8,0) 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 
             (8,2,0) 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

Running a WESTPA simulation
---------------------------

- If you encounter an issue while running the simulation

  + Use the ``--debug`` option on the servers w_run and save the output to a file.
    (note that this will generate a very detailed log of the process, try
    searching for "ERROR" for any errors and "iteration" to look at every
    iteration)
  + Use a program like hdfview, h5ls or Python with h5py library to open the
    ``west.h5`` file and ensure that the progress coordinate is being passed
    around correctly.
  + Use hdfview, h5ls or Python with h5py library to ensure that the number of
    trajectory walkers is correct.

- Is your simulation failing while the progress coordinate is being calculated?

  + One of the most error prone part during an iteration is the progress
    coordinate extraction. Programs that are not designed for quick execution
    have a lot of trouble during this step (VMD is a very commonly encountered
    one for example). Probably the best way to deal with this issue is to hard
    code a script to do the progress coordinate extraction. If you are doing
    molecular dynamics simulations multiple libraries for Python and C/C++ that
    deal with most output formats for MD packages exist and they usually come
    with a lot of convenience functions that can help you extract the progress
    coordinate. AMBER tools and GROMACS tools seems to work adequately for this
    purpose as well.

- Is your progress coordinate what you think it is?

  + Once your simulation it is running, it is well worth your time to ensure
    that the progress coordinate being reported is what you think it is. This
    can be done in a number of ways:

  + Check the ``seg_log`` output. This captures the standard error/output from
    the terminal session that your segment ran in, assuming you are running the
    executable propagator, and can be useful to ensure that everything is being
    done as you believe it should be (GROMACS tools, such as ``g_dist``, for
    instance, report what groups have their distance being calculated here).

  + Look at a structure! Do so in a program such as VMD or pyMOL, and calculate
    your progress coordinate manually and check it visually, if feasible. Does
    it look correct, and seem to match what's being reported in the .h5 file?
    This is well worth your time before the simulation has proceeded very far,
    and can save a significant amount of wallclock and computational time.

Analyzing a WESTPA simulation
-----------------------------

- If you are running the analysis on shared computing resources

  + Be sure to use the ``--serial`` flag (see the individual tool
    documentation). Otherwise, many of the included tools default to parallel
    mode (w_assign, for instance), which will create as many Python threads as
    there are CPU cores available.
