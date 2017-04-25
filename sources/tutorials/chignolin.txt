.. _chignolin_tutorial:

NAMD example: Molecular dynamics simulation of Chignolin
=========================================================

by Ernesto Suarez

updated with WESTPA version 0.9.1, changeset 3afe1b17a43f

Introduction
------------

This tutorial assumes that you are already familiar with the packages
`NAMD <http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/index.html>`__
and
`VMD <http://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/>`__
and have a basic knowledge of Python and bash scripting. Once the code
is
`installed <https://chong.chem.pitt.edu/wewiki/UserGuide:Installing>`__,
the user can
`download <https://dl.dropboxusercontent.com/u/27043352/chig_namd.tar.gz>`__
the working example directory (**chig\_namd**) where this tutorial was
run. Users have to make their own version of this folder for their
particular system. As can be noticed, there is not a unique input
file, but several key files that have to be modified. In this tutorial
we will did the modifications we needed in order to “study” the
folding in implicit solvent of a beta-hairpin, mini-protein called
`Chignolin
(GYDPETGTWG) <http://www.rcsb.org/pdb/explore/explore.do?structureId=2E4E>`__.

NAMD files needed
------------------

We will need the starting structure **seg\_initial.pdb**, and in this
case we also have the file where the velocities are stored
(**seg\_initial.vel**), since our starting point is a continuation of a
regular MD simulation. We will also need the protein structure file
**chig.psf**, the parameters file (in this case
**par\_all22\_prot\_na.prm**) and the NAMD configuration file
**md-continue.conf** that will be used to propagate the trajectories.
Note that VMD offers an automatic psf file builder from the **.pdb** via
the VMD Main menu by clicking Extensions -> Modeling -> Automatic PSF
Builder. Also notice that in any weighted ensemble simulation we do need
a stochastic component acting on the trajectories, otherwise the
trajectories would not diverge from each other. In this example we will
use Langevin dynamics (see the NAMD configuration file
**md-continue.conf** below).

md-continue.conf
~~~~~~~~~~~~~~~~

::

    # md-continue.conf
    # #   1 ps NVT production with Langevin thermostat and Generalized Born implicit solvation
    # ######################################################## INPUT #########################################################
      paraTypeCharmm                      on                      ;# Input parameter file is in 'CHARMM' format, not 'X-PLOR'
      coordinates                         parent.pdb              ;# Structure infile (.coor/.pdb)
      velocities                          parent.vel              ;# Velocity infile (.vel)
      structure                           chig.psf                ;# Topology infile (.psf)
      parameters                          par_all22_prot_na.prm   ;# Parameter infile (.prm)
    # ####################################################### ENSEMBLE #######################################################
      langevin                            on                      ;# Langevin thermostat enabled
      langevintemp                        300                     ;# Langevin thermostat temperature (K)
      langevindamping                     10                      ;# Langevin thermostat damping constant (ps -1)
      GBIS                                on                      ;# Generalized Born implicit solvent enabled
    # ################################################ NONBONDED INTERACTIONS ################################################
      cutoff                              999.0                   ;# Nonbonded cutoff for van der Waals and electrostatics (A)
      exclude                             scaled1-4               ;# Electrostatic scaling between 1-4 bonded atoms
      1-4scaling                          1.0                     ;# 1-4 interaction scaling (0.83333 for Amber force field)
      switching                           off                     ;# Van der Waals switching disabled
      pairlistdist                        1000.0                  ;# Neighbor list cut-off (A)
      stepspercycle                       250                     ;# Neighbor list update interval (timesteps)
    # ######################################################## OUTPUT ########################################################
      outputname                          seg                     ;# Output file prefix
      outputEnergies                      500                     ;# Energy log output interval (timesteps)
      restartfreq                         0                       ;# Restart file output interval (timesteps)
      binaryoutput                        no                      ;# Do not output binary trajectory
      dcdfile                             seg.dcd                 ;# Trajectory outfile (.dcd)
      dcdfreq                             500                     ;# Trajectory output interval (timesteps)
    # ###################################################### INTEGRATOR ######################################################
      timestep                            2.0                     ;# Simulation timestep (fs)
      run                                 500                     ;# Simulation durat

Other files needed
------------------

The NAMD files just mentioned are in **chig\_namd/auxfiles**. In this
folder we also have our reference file **reference.pdb** with the folded
structure, and two scripts that computes the RMSD from the reference
file, **get-rmsd-init.py** and **get-rmsd.py**. The first one is used to
initiate the simulation while the second one is used during the
propagation.(see below)

Specifying weighted ensemble parameters in the WESTPA files

The following table shows the files we will need to modify in order to
run our WESTPA simulation. In the following we will show the
modifications needed on these files in our working directory. You also
can find some general info on them
`here <https://chong.chem.pitt.edu/wewiki/UserGuide:Constructing>`__.

+------------------+--------------------------------------------------------+
| Input File       | Description                                            |
+==================+========================================================+
| env.sh           | set environment variables                              |
+------------------+--------------------------------------------------------+
| get\_pcoord.sh   | calculate progress coordinate from initial structure   |
+------------------+--------------------------------------------------------+
| system.py        | WESTPA settings (low-level)                            |
+------------------+--------------------------------------------------------+
| runseg.sh        | run segment and calculate progress coordinate          |
+------------------+--------------------------------------------------------+
| west.cfg         | WESTPA settings (high-level)                           |
+------------------+--------------------------------------------------------+
| init.sh          | initialize WESTPA files and folders                    |
+------------------+--------------------------------------------------------+
| run.sh           | run WESTPA                                             |
+------------------+--------------------------------------------------------+

env.sh
~~~~~~

This script sets environment variables that may be used across the
simulation. These include the root directory for WESTPA, the root
directory for the simulation, the name for the simulation, and the
python executable to use. It also sets the executable for NAMD; using an
environment variable for this purpose makes it easier to transition code
to different hardware or test different builds or flags of an MD code
without editing multiple files.

get\_pcoord.sh
~~~~~~~~~~~~~~

This script calculates the progress coordinate of the initial state
(**seg\_initial.pdb**), which is the RMSD from the reference structure
**reference.pdb** This script is used only during initial state
generation, where the Python script **get-rmsd-init.py** is used; during
production **runseg.sh** calculates the progress coordinate through
**get-rmsd.py**.

system.py
~~~~~~~~~

The file **system.py** contains information about the progress
coordinate, bin definitions, and the number of particles per bin. In
this example, we will be monitoring the RMSD of our protein relative to
the `folded
structure <http://www.rcsb.org/pdb/explore/explore.do?structureId=2E4E>`__,
so the number of dimensions we will consider is one:

::

    self.pcoord_ndim = 1

A system with 60 bins we can define in system.py as follows:

::

    binbounds = [0.0+0.1*i for i in xrange(60)] + [float('inf')]

We are therefore covering the entire range of our progress coordinate
since the RMSD is a finite non negative number. The bin boundaries are
left inclusive e.g. a walker with a value of 0.1 would end up in the
second bin. The positions of your bins must be either monotonically
increasing or decreasing, otherwise you will get an error message
indicating this requirement.

The number of walkers per bin is specified through the following
statement:

::

    bin.target_count = 4

We are using tau value of 1 ps, we want also to monitor the progress
coordinate every ps, writing the coordinates 1 times plus the initial
configuration:

::

    self.pcoord_len = 2

Finally, we specify the format in which the coordinates are stored:

::

    self.pcoord_dtype = numpy.float32

runseg.sh

For each segment, the WESTPA code will run the script **runseg.sh** in
each iteration. This executable creates the folders where the segment
trajectories are executed and stored, organizing them by iteration
number and segment id. In this example, runseg.sh does the following:

#. For each segment, creates a the corresponding directory
#. Copies all the files needed to run a regular NAMD molecular dynamics
   simulation into the segment directory
#. Runs dynamics using the NAMD command 'namd2' for instance:
   namd2 file.cfg > file.log
#. Computes the progress coordinate (RMSD) and records the coordinates
   to 'WEST\_COORD\_RETURN'
#. Removes any unnecessary files

west.cfg

In the configuration file called **west.cfg**, you can edit
***max\_total\_iterations*** to change the number of iterations and
***max\_run\_wallclock*** to change the maximum total time of your
simulation. If this job is being run on a computing cluster that has a
maximum wallclock time, make sure that the time specified under
**west.cfg** is less than your requested time on the cluster. Also make
sure that all directory locations are set correctly, especially under
the ***data\_refs*** heading in the ***data*** section, and throughout
the ***executable*** section since this is how WE code knows where to
send the output, where the initial states are, where the progress
coordinate lies, etc. Everything should already be set correctly, but
make sure that the variable 'WEST\_SIM\_ROOT', located in the **env.sh**
file, is set to your where your simulation directory exists.

init.sh

This script initializes the WESTPA system. It removes files from
previous runs, and uses get\_pcoord.sh to generate initial states

run.sh

This script is used to run WESTPA.

Running the simulation

From the WEST\_SIM\_ROOT (simulation root) directory, you can initialize
the simulation by entering the following at the command prompt:

::

    ./init.sh

The script should create three directories, **traj\_segs**, and
**seg\_logs**, as well as an HDF5 file named **west.h5**. It calls
**w\_init.py** from the main weighted ensemble code.

Now that your simulation has been initialized, it is ready to be run by
the weighted ensemble code. Use the command:

::

    ./run.sh

The script **run.sh** calls **w\_run.py** from the main weighted
ensemble code. If it doesn't work, check to see if the **env.sh** is set
up properly and if it points to the right directory for your weighted
ensemble code. Make sure that the 'WEST\_ROOT' variable is set to where
the **westpa** directory exists and the 'WEST\_SIM\_ROOT' variable is
set to where your simulation directory exists.

If this simulation is being run on a computing cluster, **w\_run** can
be executed from a batch script. See the clusters page for more
information on how to submit jobs to specific clusters.

Restarting the simulation
-------------------------

In **west.cfg** we have set ***max\_total\_iterations*** equals to 300.
We might increase this number up to the desired value and run our
simulation again. The westpa code will continue the simulation from
where you left off, based on the data present in the **west.h5** file
(see below). If you wanted to restart the simulation from the beginning,
you would need to run the **init.sh** script again, which would remove
the existing **west.h5** file and create a new one. Once you have
changed the ***max\_total\_iterations*** flag to, say 400, execute the
**run.sh** script again.

Analyzing the data
------------------

Output
------

The results of the simulations are stored in the **traj\_segs** folder.
Here you can find the results of each small segment. It is organized by
iterations, within which are directories for each segment (as specified
in **west.cfg**). And within each segment directory are all of the NAMD
output files. Because we specified 300 iterations in the **west.cfg**
file, you should see 300 iteration directories, named **000001**,
**000002**...

In **west.h5** is stored the most relevant information about our
simulation. It is convenient for the user to download the program
`HDFView <http://www.hdfgroup.org/hdf-java-html/hdfview/>`__, that allow
us to see how the information is stored and helps us to create tools for
the analysis.

Observing folding events

Through the following python script, we can monitor the progress of the
simulation by plotting the C-alpha RMSD (in A) of the
closest/farest segment to the folded structure vs. the total simulation
time.

::

    #!/usr/bin/env python

    import sys
    import numpy as np
    import h5py

    if(len(sys.argv)<2):
          print "Usage:"
          print "getMinMax.py <argv1>"
          print "<argv1>: Num of Iterations"
          sys.exit(-1)

    totNumIter=int(sys.argv[1])

    try:
        datafile=h5py.File("west.h5","r")
    except IOError as err:
        print err.errno
        print err.strerror
        sys.exit(-1)

    counter = 0
    for I in xrange(1,totNumIter+1):
          itr = "iter_"+str(I).rjust(8,'0')
          coords = np.array(datafile["iterations"][itr]["pcoord"])
          numOfsegments = len(coords)
          counter += numOfsegments      #Stores the total sim time in Tau units
          firstcoord = coords[:,1,0]
          firstcoord.sort()
          print counter, firstcoord[0],firstcoord[-1]

In this example, we define the folded state as any structure with a
C-alpha RMSD less than 2.5A. As shown in the plot, folding events
occur after the first 20ns of simulation.

.. figure:: ../_static/chig_minmax.png

Computing the folding rate constant

WESTPA includes several scripts for analysis located in
**$WEST\_ROOT/bin**. In *init.sh* we have specified the folded structure
as those structures with RMSD less than 2.5A. Using
**w\_fluxanl**, we can calculate the flux into this target (folded)
state, and from that calculate the folding rate. **w\_fluxanl** may be
run with the following commands:

::

    source env.sh
    $WEST_ROOT/bin/w_fluxanl

The script will output the flux into the target state including
confidence intervals calculated using the block bootstrap method.

::

    Calculating mean flux and confidence intervals for iterations [1, 301)
    target 'bound':
      correlation length = 0 tau
      mean flux and CI   = 5.028992e-03 (4.231789e-04,1.185007e-02) tau^(-1)

Taking the inverse of this flux (1 / 5.028992e-03) yields the folding
rate in units of τ, the segment duration.

More information on how to use *w\_fluxanl* can be viewed using the
'--help' flag. *w\_fluxanl* also stores this information in an hdf5
file, *fluxanl.h5*. Using the python libraries h5py and pylab, we can
visualize this data. Open a python interpreter and run the following
commands:

::

    import h5py, numpy, pylab
    fluxanl              = h5py.File('fluxanl.h5')
    flux                 = numpy.zeros(300)
    first_binding        = 300 - fluxanl['target_flux']['target_0']['flux'].shape[0]
    flux[first_binding:] = numpy.array(fluxanl['target_flux']['target_0']['flux'])
    pylab.xlim([289,300])
    pylab.plot(flux)
    pylab.xlabel("Iteration")
    pylab.ylabel("Instantaneous Flux $(\\frac{1}{\\tau})$")
    pylab.show()

.. figure:: ../_static/chig_flux1.png

The x-axis represents the iteration number, and the y-axis the flux into
the bound state in units of τ\ :sup:`-1`. For my simulation, the first
binding event occurred in iteration 289. The instantaneous flux is noisy
and can be difficult to interpret, so let's plot the time evolution of
flux as well. Run **w\_fluxanl** again, this time with the '--evol'
flag.

::

    $WEST_ROOT/bin/w_fluxanl --evol

This will add a dataset named ['flux\_evolution'] to the ['target\_0']
group in 'fluxanl.h5'. We may plot the time evolution of flux using the
following commands at a python interpreter:

::

    import h5py, numpy, pylab
    fluxanl   = h5py.File('fluxanl.h5')
    mean_flux = numpy.zeros(300)
    ci_ub     = numpy.zeros(300)
    ci_lb     = numpy.zeros(300)
    first_binding             = 300 - fluxanl['target_flux']['target_0']['flux_evolution']['mean'].shape[0]
    mean_flux[first_binding:] = numpy.array(fluxanl['target_flux']['target_0']['flux_evolution']['mean'])
    ci_lb[first_binding:]     = numpy.array(fluxanl['target_flux']['target_0']['flux_evolution']['ci_lb'])
    ci_ub[first_binding:]     = numpy.array(fluxanl['target_flux']['target_0']['flux_evolution']['ci_ub'])
    pylab.plot(mean_flux, 'b', ci_lb, 'g', ci_ub, 'r')
    pylab.xlim([289,300])
    pylab.xlabel("Iteration")
    pylab.ylabel("Mean Flux $(\\frac{1}{\\tau})$")
    pylab.show()

.. figure:: ../_static/chig_flux2.png

Since our simulation is relatively short, we do not expect to observe a
constant flux and the end of the simulation. The mean flux obtained from
the output of **w\_fluxanl** for our simulation was 5.03 x 10\ :sup:`-3`
in units of τ\ :sup:`-1`. Since τ for our simulation was 1ps, the
folding rate is 5.03 x 10\ :sup:`-3` ps\ :sup:`-1` with a 95% CI of 4.2
x10\ :sup:`-4` to 1.18 x10\ :sup:`-2`. In order to obtain a more precise
association rate and a relatively constant flux, the user would need to
run more iterations of the simulation, which may easily be done by
changing *west.cfg*.

Useful links
------------

-  `NAMD tutorial from the official web
   page <http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/index.html>`__,
-  `VMD tutorial from the official web
   page <http://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/>`__
