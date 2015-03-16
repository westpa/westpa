OpenMM Tutorial: Molecular Dynamics of Na\ :sup:`+`/Cl\ :sup:`-` Association
============================================================================

by Joshua L. Adelman

Updated with WESTPA version 1.0 beta and OpenMM 5.1

Overview
--------

Requirements: ~? hours wallclock time on ?; ~? GB disk space

In this tutorial we will use the standard weighted ensemble approach to
simulate Na\ :sup:`+`/Cl\ :sup:`-` association in Generalized Born implicit
solvent. The system consists of single Na\ :sup:`+` and Cl\ :sup:`-` ions
modeled with the
`Amber force field <http://ambermd.org/#ff>`_,
using the distance between the two ions as the progress coordinate.
`OpenMM <http://openmm.org>`_ will be used to run the
molecular dynamics, and familiarity with it is a prerequisite (see `tutorials
<http://docs.openmm.org/6.2.0/userguide/index.html>`_).

This tutorial uses the same starting files, generated using Amber, as the
:ref:`Introductory Amber Tutorial <amber_tutorial>`.
Instead of using Amber (``pmemd`` or ``sander``) to run the dynamics, we
will use the OpenMM framework along with its python wrapper to propagate the
system.

While the Amber tutorial uses the executable propagator to call out to Amber
and then a number of shell scripts to communicate the results back to WESTPA,
OpenMM's python wrapper allows us to integrate the molecular dynamics
simulation directly into WESTPA. This tutorial will therefore require some
understanding of python and the OpenMM framework. Additionally, the focus will
be on integrating OpenMM into WESTPA and less on the analysis of the resulting
simulations.

Since the Na\ :sup:`+`/Cl\ :sup:`-` system only contains two atoms, we will
make some design decisions in setting up the simulation that may not be
appropriate more generally. Most notably, we will avoid writing trajectory data
and restart files for each trajectory segment to disk as separate files, as is
generally done when using the executable propagator. Instead, we will store all
of this data in the main hdf5 file that also contains the data related to the
weighted ensemble run.

The first step is to set up a directory containing the necessary AMBER and
WESTPA files. A working example directory can be found at
``westpa/lib/examples/nacl_openmm``.

Preparing the OpenMM files
--------------------------

======================= ===============================================
Input File              Description
======================= ===============================================
...                     ...
======================= ===============================================

We will begin with the coordinate and topology files (``nacl.inpcrd`` and
``nacl.prmtop`` respectively) produced by AmberTools. We will construct an
OpenMM system and integrator and then serialize the resulting python objects
using a short script::

  python build_system.py

Preparing the WESTPA files
--------------------------

======================= =======================================================
Input File              Description
======================= =======================================================
env.sh                  set environment variables
system.py               system implementation
openmm_propagator.py    custom propagator to run openmm calculation across
                        multiple devices
restart_plugin.py       plugin to allow restart information to be stored in the
                        west.h5 file
west.cfg                WESTPA configuration
init.sh                 initialize WESTPA
======================= =======================================================

system.py
~~~~~~~~~

This file contains information about the progress coodinate, binning, walkers
per bin and more. This file is nearly identical to one defined in the
`Introductory Amber Tutorial <amber_tutorial>`. In this example
we will be using the distance between the two ions as the progress coordinate,
giving us a one dimensional coordinate::

  self.pcoord_ndim = 1

The positions of the bins along this progress coordinate are the same as those
used by `Zwier, Kaus, and Chong
<http://pubs.acs.org/doi/abs/10.1021/ct100626x>`_::

  binbounds = [0.0] + [2.8, 2.88, 3.0, 3.10, 3.29, 3.79, 3.94, 4.12, 4.39,
    5.43] + [5.90+1.0*i for i in xrange(0,11)] + [30,float('inf')]

Since every walker must lie in a bin, the upper boundary to the last bin is set
to infinity i.e. ``[30, float('inf')]``. The bin boundaries are left inclusive
e.g. a walker with a value of 2.8 would end up in the second bin. The positions
of your bins must be either monotonically increasing or decreasing - otherwise,
you will get an error message indicating this requirement.

The number of walkers per bin is specified by the following::

  bin.target_count = 48

Using a tau value of 0.5 ps, we will monitor the progress coordinate every 0.05
ps, writing coordinates 10 times. Including the initial configuration this
gives an expected progress coordinate length of 11:

  self.pcoord_len = 11

Finally, we specify the format in which the coordinates are stored:

  self.pcoord_dtype = numpy.float32

openmm_propagator.py
~~~~~~~~~~~~~~~~~~~~~

The OpenMMPropagator subclasses the WESTPropagator interface and implements all
of the necessary methods to run a WESTPA simulation using OpenMM. The
implementation presented in this example, while fairly generic, is still
specific enough to the Na\ :sup:`+`/Cl\ :sup:`-` association example, that
changes will likely be necessary to adapt it for another system. Below is a
brief description of each method in the class:

__init__ method
_______________

The ``__init__`` method is primarily responsible for parsing the configuration
parameters form ``west.cfg`` and building the OpenMM system, integrator and
platform objects. Since each OpenMM context must be tied to a unique
integrator, the ``propagator`` method actually deserializes the integrator for
each propagation step. In this method, however, it is primarily being used to
retrieve the temperature of the system.

static methods
______________

The OpenMMPropagator contains three methods that are tagged with the Python's
``@staticmethod`` decorator. This designation just allows the methods to be
encapsulated within the class, but they do not have direct access to the
class's internal data. The ``dist`` method just calculates a simple Euclidean
distance between two points and is used in calculating the pcoord of a
conformation of the system. The ``makepath`` method assembles a path on the
filesystem from a template and is used to tell the propagator where to grab
initial state information from. The ``mkdir_p`` method augments the standard
library's ``os`` module to allow unix ``mkdir -p`` like behavior.

get_pcoord method
_________________

This method assigns a pcoord value to a given state. The state can either be an
``BasisState``, in which case we uses the basis state's coordinate, which are
stored as a class variable to calculate the pcoord using ``dist``. If the state
is an ``InitialState`` (i.e. the result of perturbing the x-position of one of
the ions by a random amount), we construct the path to the file containing its
coordinates, and calculate the pcoord after reading the file from disk.

propagate method
________________

The ``propagate`` method takes a set of segments and runs each for a length of
time tau. Initially, the method attempts to assign the calculation to a device
based on the ``WM_PROCESS_INDEX`` environment variable if it is available (both
the zmq and processes work managers set it, but the other work managers do
not). A context is then constructed, before the method iterates over all
segments.

For each segment, an initial set of coordinates or velocities are obtained
either from the parent segment, if this segment is a continuation of previous
dynamics, or from an initial state if the segment is being initiated at the
start of the WE calculation or is the result of a recycling event. Dynamics are
then run using the OpenMM integrator. At a user-specified interval, the
calculation is halted and the coordinates and velocities, along with the
calculated pcoord are saved to temporary arrays. Finally this data is
transferred to the segment's internal data structures.

gen_istate method
_________________

This method takes a basis state and generates an initial state by randomly
perturbing the basis state and storing the results to disk using the naming
convention specified by the template given in the ``west.cfg`` file.

restart_plugin.py
~~~~~~~~~~~~~~~~~~

In order to restart a segment from its parent, we need access to the last set
of coordinates and velocities recorded for the parent in the ``coord`` and
``veloc`` data sets. We use a custom plugin that is run just before the
propagation step that temporarily loads the necessary coordinates and
velocities into a segment's data dictionary as
``segment.data['restart_coord']`` and ``segment.data['restart_veloc']``. The
propagator will then delete this data once it has been transferred to the
OpenMM context.

This allows us to run the entire simulation from the main hdf5 file without
writing any per-segment data to individual files. While convenient for a simple
system like the one in this example, it may not be as desirable for systems
with a large number of particles. In that case the propagator will need to be
modified to load the restart data from individual files contained in the
traj_segs directory on the file system, as is the case for the examples that
use the executable propagator.

west.cfg
~~~~~~~~

The actual WESTPA simulation is configured using the yaml-formatted
``west.cfg`` file. The custom propagator will extract a number of parameters
from the ``openmm`` section shown below.::

  ---
  west:
    ...
    openmm:
      system: 
        file: system.xml
      integrator: 
        file: integrator.xml
        steps_per_tau: 250
        steps_per_write: 25
      platform:
        name: CUDA
        #properties: {'OpenCLPlatformIndex': '1', 'OpenCLDeviceIndex': '0'} # Platform specific properties 

The xml files are the output of running the ``build_system.py`` script. Within
the ``integrator`` section, the ``steps_per_tau`` and ``steps_per_write``
specify the number of time steps that the integrator should advance the system
per tau (so 250 x 2 fs = 0.5 ps) and at what frequency, in numbers of steps,
that the pcoord and auxiliary data should be collected, respectively.

The ``platform`` section defines a platform ``name``, which can be
``Reference``, ``CUDA``, or ``OpenCL``, assuming the latter two are installed
on your system. The CUDA platform requires a compatible GPU card, but the
OpenCL platform, in addition to running on GPUs supports both the Intel and AMD
CPU OpenCL SDK.

Finally, the ``properties`` variable under the ``platform`` section defines a
dictionary, whose members override the defaults specified in the propagator
``__init__`` method. See the defaults for all possible platform specific
settings. Importantly, the ``XXXDeviceIndex`` settings are ignored when running
in parallel using either the zeromq or processes work managers, since they set
that variable dynamically for each worker. However, when running in serial mode
on a multi-device system, it can be useful to select a specific device to run
the calculation on. When running using the OpenCL platform, the `oclutils
<https://github.com/nbigaouette/oclutils>`_ library is useful in extracting
information about the available devices and platforms (in the OpenCL meaning of
platform, rather than the OpenMM one).

There are also some important settings under the ``propagation`` section::

  ---
  west:
    ...
    propagation:
      max_total_iterations: 2
      max_run_wallclock: 2:00:00
      propagator: openmm_propagator.OpenMMPropagator
      gen_istates: true
      block_size: 138

In addition to setting the location of the custom openmm propagator, this
section allows you to set the total number of iterations to run using
``max_total_iterations``. This should be changed to collect data for this
system to at least 100. The ``max_run_wallclock`` time should also be adjusted
depending on the hardware being used to run this simulation. Using four GTX
680s, this system takes approximately 16 seconds per iteration.

A particularly important setting in terms of the performance of the calculation
is ``block_size``. This parameter determines how many segments are sent to the
propagator at a time during the run. Since setting up the OpenMM context is
quite expensive, one can get a large boost in performance by re-using the same
context and just pushing new coordinates and velocities to it. So if the
calculation is run using the serial work manager, ``block_size`` should be set
to the maximum number of replicas possible for the system, which in this case
is 552. Likewise, if running the calculation over 4 devices, this number should
be 552 / 4 = 138.

Running the simulation
----------------------

The simulation can then be initiated and ran using the shell scripts,
``init.sh`` and ``run.sh``.

From the simulation root directory (``$WEST_SIM_ROOT``) directory, enter into
the command line::

    ./init.sh

The script should create a directory called ``istates``, as well as an HDF5
file named ``west.h5``. Because the ``gen_istates`` flag was set to True in the
west.cfg file, the propagator's ``gen_istate`` method should prepare multiple
different ``.txt`` input coordinate files, located in the ``istates/1``
directory. The ``init.sh`` script should finish by printing "Simulation
prepared." with a short list (8 lines) of probabilities and statistics about
the initial state of the methane-methane simulation.

Now that your simulation has been initialized, it is ready to be run by the
weighted ensemble code. Use the command::

  ./run.sh --work-manager=zmq --n-workers=4 &

to use the ``zmq`` work manager and run using 4 workers.

The ``init.sh`` and ``run.sh`` scripts call ``w_init.py`` and ``w_run.py`` from
the main weighted ensemble code, respectively. If either does not work, check
to see if the ``env.sh`` is set up properly and if it points to the right
directory for your weighted ensemble code (the default settings assume you are
running from within the westpa/lib/nacl_openmm directory). Make sure that the
``WEST_ROOT`` variable is set to where the ``westpa`` directory exists and the
``WEST_SIM_ROOT`` variable is set to where your simulation directory exists.

Analyzing the data
------------------

Output
~~~~~~

======================= =======================================================
Output File             Remarks
======================= =======================================================
west.h5                 WESTPA output in hdf5 database
west.log                WESTPA log file
======================= =======================================================

The way in which we set up the calculation, all output data is stored within
the hdf5 file, ``west.h5``. Because we specified 2 iterations in the
``west.cfg`` file, the simulation should have only run for a short period of
time. This is not enough to generate any meaningful results, but is sufficient
to ensure that the system was set up properly.

In the ``west.cfg`` file, change the ``max_total_iterations`` variable to 100.
The westpa code will continue the simulation from where you left off, based on
the data present in the ``west.h5`` file. If you wanted to restart the
simulation from scratch, you would need to run the ``init.sh`` script again,
which would remove the existing ``west.h5`` file and create a new one. Once you
have changed the ``max_total_iterations`` flag to 100, execute the ``run.sh``
script again. Simulating 100 iterations may take some time, so be prepared to
wait. Using 4 GTX 680s and running with the CUDA, platform, this should take
about 25 minutes. Not, that for a small number of atoms, such is the case for
this system, running on the GPUs does not leverage the full capabilities of the
hardware and is likely to be slower than using an optimized CPU-based code.

Computing the association rate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WESTPA includes several tools for analysis located in ``$WEST_ROOT/bin``. In
``init.sh`` we specified the bin containing an Na\ :sup:`+`/Cl\ :sup:`-`
distance of 1.8 Å as the bound state, and that containing a distance of 16.9 Å
as the unbound state. Using ``w_fluxanl``, we can calculate the flux into these
target states, and from that calculate the association rate of Na\ :sup:`+`/Cl\
:sup:`-`. ``w_fluxanl`` may be run with the following commands::

  source env.sh
  $WEST_ROOT/bin/w_fluxanl

The script will output the flux into the target states including confidence
intervals calculated using the block bootstrap method::

  Calculating mean flux and confidence intervals for iterations [1,101)
  target 'bound':
    correlation length = w tau
    mean flux and CI   = x (y, z) tau^(-1)

More information on how to use ``w_fluxanl`` can be viewed using the ``--help``
flag. ``w_fluxanl`` also stores this information in an hdf5 file,
``fluxanl.h5``.

Presently, ``w_fluxanl`` has used the data from all 100 iterations (note the
exclusive bracket after 101) to calculate the mean flux (x) and the 95%
confidence interval (y, z) for reaching the bound state (target 'bound'), which
we specified as less than 2.8 angstroms of separation in the ``system.py`` file
and with the target state variable in ``init.sh``. The value given for the flux
also represents the association rate. Taking the inverse of the mean flux (1/x)
will give the mean first passage time for Na\ :sup:`+`/Cl\ :sup:`-` in units of
tau. We can further analyze the output of ``w_fluxanl`` by investigating the
``fluxanl.h5`` file. You can look at the data contained within the file by
using programs such as h5ls or hdfview, but I am instead going to use h5py in
python to analyze the data. Open up ``ipython`` in the interactive plotting
mode::

  ipython --pylab

and then enter the following commands::

  import h5py
  import numpy as np

  fluxanl = h5py.File('fluxanl.h5')
  fluxanl['target_flux']['index'][:]

We can see that the dataset named ['index'] contains the output printed
above by ``w_fluxanl``. We can plot the flux using::

  flux = np.array(fluxanl['target_flux']['target_0']['flux'])
  plot(flux)

The x-axis represents the iteration number recorded after the occurence of the
first binding event. The y-axis represents the flux in units of tau\ :sup:`-1`.
We can see that the instantaneous flux has settled after large fluctuations
during the first part of the run, however the plot is also relatively noisy. To
reduce noise, we can plot the time evolution flux. Run the ``w_fluxanl`` tool
again, this time with the '--evol' flag at the end of the command. Running this
command will add an HDF5 dataset named ['flux_evolution'] to the ['target_0']
group. To plot the time evolution flux, you can use the following python code,
continuing from the above ipython session::

  mean_flux = fluxanl['target_flux']['target_0']['flux_evolution']['expected']
  ci_lb = fluxanl['target_flux']['target_0']['flux_evolution']['ci_lbound']
  ci_ub = fluxanl['target_flux']['target_0']['flux_evolution']['ci_ubound']
  plot(mean_flux, 'b', ci_lb, 'g', ci_ub, 'r')

Compared to the first plot of the instantaneous flux, the time evolution
plot is much less noisy. We can see that the flux is leveling off and
the confience intervals have somewhat converged, meaning that the
simulation is approaching steady-state conditions.

Visualizing a selected pathway
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to visualize a binding event, you will need to stitch together the
individual trajectory segments that start at the initial state and then reach
the bound state. The `introductory Amber tutorial <amber_tutorial>` provides
directions on how to extract the sequence of segments in a set of successful
binding events, however the script to construct a visualization of the pathway
will not work for this example since we have stored all of the relevant data
directly in the ``west.h5`` file. For this example, we leave writing the
necessary script as an exercise. To create a netcdf-formatted Amber trajectory
file, you might want to take a look at `netcdf4storage.py
<https://bitbucket.org/joshua.adelman/stringmethodexamples/raw/tip/shared/elasticnetwork-langevin/netcdf4storage.py>`_
or you might consider using the dcd writer built into OpenMM which can imported
into python using::

  import simtk.openmm.app.dcdfile


Useful links
------------

- `Official OpenMM web page <http://openmm.org>`_
- `OpenMM tutorials from the official web page
  <http://docs.openmm.org/6.2.0/userguide/index.html>`_

Useful hints
------------

- Make sure your paths are set correctly in ``env.sh``
- If the simulation doesn't stop properly with CTRL+C , use CTRL+Z.
- Another method to stop the simulation relatively cleanly is to rename
  ``runseg.sh``; WESTPA will shut the simulation down and prevent the hdf5 file
  from becoming corrupted. Some extra steps may be necessary to ensure that the
  analysis scripts can be run successfully.

References
----------

- `Zwier, MC, Kaus, JW, Chong, LT. Efficient Explicit-Solvent Molecular
  Dynamics Simulations of Molecular Association Kinetics: Methane/Methane,
  Na+/Cl−, Methane/Benzene, and K+/18-Crown-6 Ether. J Chem Theory Comput.
  2011. <http://pubs.acs.org/doi/abs/10.1021/ct100626x>`_
