.. _setup:

Setup
=====

Defining and Calculating Progress Coordinates
---------------------------------------------

Binning
-------

The Weighted Ensemble method enhances sampling by partitioning the space
defined by the progress coordinates into non-overlapping bins. WESTPA provides
a number of pre-defined types of bins that the user must parameterize within
the system.py file, which are detailed below.

Users are also free to implement their own mappers. A bin mapper must
implement, at least, an ``assign(coords, mask=None, output=None)`` method,
which is responsible for mapping each of the vector of coordinate tuples
``coords`` to an integer (``numpy.uint16``) indicating what bin that coordinate
tuple falls into. The optional ``mask`` (a numpy bool array) specifies that
some coordinates are to be skipped; this is used, for instance, by the
recursive (nested) bin mapper to minimize the number of calculations required
to definitively assign a coordinate tuple to a bin. Similarly, the optional
``output`` must be an integer (``uint16``) array of the same length as
``coords``, into which assignments are written. The ``assign()`` function must
return a reference to ``output``. (This is used to avoid allocating many
temporary output arrays in complex binning scenarios.)

A user-defined bin mapper must also make an ``nbins`` property available,
containing the total number of bins within the mapper.

RectilinearBinMapper
~~~~~~~~~~~~~~~~~~~~

Creates an N-dimensional grid of bins. The Rectilinear bin mapper is
initialized by defining a set of bin boundaries:

::

  self.bin_mapper = RectilinearBinMapper(boundaries)

where ``boundaries`` is a list or other iterable containing the bin boundaries
along each dimension. The bin boundaries must be monotonically increasing along
each dimension. It is important to note that a one-dimensional bin space must
still be represented as a list of lists as in the following example:::

  bounds = [-float('inf'), 0.0, 1.0, 2.0, 3.0, float('inf')]
  self.bin_mapper = RectilinearBinMapper([bounds])

A two-dimensional system might look like:::

  boundaries = [(-1,-0.5,0,0.5,1), (-1,-0.5,0,0.5,1)] 
  self.bin_mapper = RectilinearBinMapper(boundaries)

where the first tuple in the list defines the boundaries along the first
progress coordinate, and the second tuple defines the boundaries along the
second. Of course a list of arbitrary dimensions can be defined to create an
N-dimensional grid discretizing the progress coordinate space.

VoronoiBinMapper
~~~~~~~~~~~~~~~~

A one-dimensional mapper which assigns a multidimensional progress coordinate
to the closest center based on a distance metric. The Voronoi bin mapper is
initialized with the following signature within the
``WESTSystem.initialize``:::

  self.bin_mapper = VoronoiBinMapper(dfunc, centers, dfargs=None, dfkwargs=None)

- ``centers`` is a ``(n_centers, pcoord_ndim)`` shaped numpy array defining
  the generators of the Voronoi cells
- ``dfunc`` is a method written in Python that returns an ``(n_centers, )``
  shaped array containing the distance between a single set of progress
  coordinates for a segment and all of the centers defining the Voronoi
  tessellation. It takes the general form:::

    def dfunc(p, centers, *dfargs, **dfkwargs):
        ...
        return d

where ``p`` is the progress coordinates of a single segment at one time slice
of shape ``(pcoord_ndim,)``, ``centers`` is the full set of centers, ``dfargs``
is a tuple or list of positional arguments and ``dfwargs`` is a dictionary of
keyword arguments. The bin mapper's ``assign`` method then assigns the progress
coordinates to the closest bin (minimum distance). It is the responsibility of
the user to ensure that the distance is calculated using the appropriate
metric.

- ``dfargs`` is an optional list or tuple of positional arguments to pass into
  ``dfunc``.
- ``dfkwargs`` is an optional dict of keyword arguments to pass into ``dfunc``.

FuncBinMapper
~~~~~~~~~~~~~

A bin mapper that employs a set of user-defined function, which directly
calculate bin assignments for a number of coordinate values. The function is
responsible for iterating over the entire coordinate set. This is best used
with C/Cython/Numba methods, or intellegently-tuned numpy-based Python
functions.

The ``FuncBinMapper`` is initialized as:::

  self.bin_mapper = FuncBinMapper(func, nbins, args=None, kwargs=None)

where ``func`` is the user-defined method to assign coordinates to bins,
``nbins`` is the number of bins in the partitioning space, and ``args`` and
``kwargs`` are optional positional and keyword arguments, respectively, that
are passed into ``func`` when it is called.

The user-defined function should have the following form:::

  def func(coords, mask, output, *args, **kwargs)
      ....

where the assignments returned in the ``output`` array, which is modified
in-place.

As a contrived example, the following function would assign all segments to bin
0 if the sum of the first two progress coordinates was less than ``s*0.5``, and
to bin 1 otherwise, where ``s=1.5``:::

  def func(coords, mask, output, s):
      output[coords[:,0] + coords[:,1] < s*0.5] = 0
      output[coords[:,0] + coords[:,1] >= s*0.5] = 1
   
  ....

  self.bin_mapper = FuncBinMapper(func, 2, args=(1.5,)) 

VectorizingFuncBinMapper
~~~~~~~~~~~~~~~~~~~~~~~~

Like the ``FuncBinMapper``, the ``VectorizingFuncBinMapper`` uses a
user-defined method to calculate bin assignments. They differ, however, in that
while the user-defined method passed to an instance of the ``FuncBinMapper`` is
responsible for iterating over all coordinate sets passed to it, the function
associated with the ``VectorizingFuncBinMapper`` is evaluated once for each
unmasked coordinate tuple provided. It is not responsible explicitly for
iterating over multiple progress coordinate sets.

The ``VectorizingFuncBinMapper`` is initialized as:::

  self.bin_mapper = VectorizingFuncBinMapper(func, nbins, args=None, kwargs=None)

where ``func`` is the user-defined method to assign coordinates to bins,
``nbins`` is the number of bins in the partitioning space, and ``args`` and
``kwargs`` are optional positional and keyword arguments, respectively, that
are passed into ``func`` when it is called.

The user-defined function should have the following form:::

  def func(coords, *args, **kwargs)
      ....

Mirroring the simple example shown for the ``FuncBinMapper``, the following
should result in the same result for a given set of coordinates. Here segments
would be assigned to bin 0 if the sum of the first two progress coordinates was
less than ``s*0.5``, and to bin 1 otherwise, where ``s=1.5``:::

  def func(coords, s):
      if coords[0] + coords[1] < s*0.5:
          return 0
      else:
          return 1
  ....

  self.bin_mapper = VectorizingFuncBinMapper(func, 2, args=(1.5,)) 

PiecewiseBinMapper
~~~~~~~~~~~~~~~~~~

RecursiveBinMapper
~~~~~~~~~~~~~~~~~~

The ``RecursiveBinMapper`` is used for assembling more complex bin spaces from
simpler components and nesting one set of bins within another. It is
initialized as:::

  self.bin_mapper = RecursiveBinMapper(base_mapper, start_index=0)

The ``base_mapper`` is an instance of one of the other bin mappers, and
``start_index`` is an (optional) offset for indexing the bins. Starting with
the ``base_mapper``, additional bins can be nested into it using the
``add_mapper(mapper, replaces_bin_at)``. This method will replace the bin
containing the coordinate tuple ``replaces_bin_at`` with the mapper specified
by ``mapper``.

As a simple example consider a bin space in which the ``base_mapper`` assigns a
segment with progress coordinate with values <1 into one bin and >= 1 into
another. Within the former bin, we will nest a second mapper which partitions
progress coordinate space into one bin for progress coordinate values <0.5 and
another for progress coordinates with values >=0.5. The bin space would look
like the following with corresponding code:::

  '''         
               0                            1                      2
               +----------------------------+----------------------+
               |            0.5             |                      |
               | +-----------+------------+ |                      |
               | |           |            | |                      |
               | |     1     |     2      | |          0           |
               | |           |            | |                      |
               | |           |            | |                      |
               | +-----------+------------+ |                      |prettyprint
               +---------------------------------------------------+      
  '''

  def fn1(coords, mask, output):
      test = coords[:,0] < 1
      output[mask & test] = 0
      output[mask & ~test] = 1
    
  def fn2(coords, mask, output):
      test = coords[:,0] < 0.5
      output[mask & test] = 0
      output[mask & ~test] = 1

  outer_mapper = FuncBinMapper(fn1,2)
  inner_mapper = FuncBinMapper(fn2,2)
  rmapper = RecursiveBinMapper(outer_mapper)
  rmapper.add_mapper(inner_mapper, [0.5])

Examples of more complicated nesting schemes can be found in the `tests
<https://github.com/westpa/westpa/blob/master/lib/west_tools/tests/testbinning.py>`_
for the WESTPA binning apparatus.

Initial/Basis States
--------------------

A WESTPA simulation is initialized using ``w_init`` with an initial
distribution of replicas generated from a set of basis states. These basis
states are used to generate initial states for new trajectories, either at the
beginning of the simulation or due to recycling. Basis states are specified
when running ``w_init`` either in a file specified with ``--bstates-from``, or
by one or more ``--bstate`` arguments. If neither ``--bstates-from`` nor at
least one ``--bstate`` argument is provided, then a default basis state of
probability one identified by the state ID zero and label "basis" will be
created (a warning will be printed in this case, to remind you of this
behavior, in case it is not what you wanted).

When using a file passed to ``w_init`` using ``--bstates-from``, each line in
that file defines a state, and contains a label, the probability, and
optionally a data reference, separated by whitespace, as in:::

  unbound    1.0
  
or::

  unbound_0    0.6        state0.pdb
  unbound_1    0.4        state1.pdb

Basis states can also be supplied at the command line using one or more
``--bstate`` flags, where the argument matches the format used in the state
file above. The total probability summed over all basis states should equal
unity, however WESTPA will renormalize the distribution if this condition is
not met.

Initial states are the generated from the basis states by optionally applying
some perturbation or modification to the basis state. For example if WESTPA was
being used to simulate ligand binding, one might want to have a basis state
where the ligand was some set distance from the binding partner, and initial
states are generated by randomly orienting the ligand at that distance. When
using the executable propagator, this is done using the script specified under
the ``gen_istate`` section of the ``executable`` configuration. Otherwise, if
defining a custom propagator, the user must override the ``gen_istate`` method
of ``WESTPropagator``.

When using the executable propagator, the the script specified by
``gen_istate`` should take the data supplied by the environmental variable
``$WEST_BSTATE_DATA_REF`` and return the generated initial state to
``$WEST_ISTATE_DATA_REF``. If no transform need be performed, the user may
simply copy the data directly without modification. This data will then be
available via ``$WEST_PARENT_DATA_REF`` if ``$WEST_CURRENT_SEG_INITPOINT_TYPE``
is ``SEG_INITPOINT_NEWTRAJ``.

Target States
-------------

WESTPA can be run in a recycling mode in which replicas reaching a target state
are removed from the simulation and their weights are assigned to new replicas
created from one of the initial states. This mode creates a non-equilibrium
steady-state that isolates members of the trajectory ensemble originating in
the set of initial states and transitioning to the target states. The flux of
probability into the target state is then inversely proportional to the mean
first passage time (MFPT) of the transition.

Target states are defined when initializing a WESTPA simulation when calling
``w_init``. Target states are specified either in a file specified with
``--tstates-from``, or by one or more ``--tstate`` arguments. If neither
``--tstates-from`` nor at least one ``--tstate`` argument is provided, then an
equilibrium simulation (without any sinks) will be performed.

Target states can be defined using a text file, where each line defines a
state, and contains a label followed by a representative progress coordinate
value, separated by whitespace, as in:::

  bound     0.02

for a single target and one-dimensional progress coordinates or:::

  bound    2.7    0.0
  drift    100    50.0

for two targets and a two-dimensional progress coordinate.

The argument associated with ``--tstate`` is a string of the form ``'label,
pcoord0 [,pcoord1[,...]]'``, similar to a line in the example target state
definition file above. This argument may be specified more than once, in which
case the given states are appended to the list of target states for the
simulation in the order they appear on the command line, after those that are
specified by ``--tstates-from``, if any.

WESTPA uses the representative progress coordinate of a target-state and
converts the **entire** bin containing that progress coordinate into a
recycling sink.

Propagators
-----------

The Executable Propagator
~~~~~~~~~~~~~~~~~~~~~~~~~

Writing custom propagators
~~~~~~~~~~~~~~~~~~~~~~~~~~

While most users will use the Executable propagator to run dynamics by calling
out to an external piece of software, it is possible to write custom
propagators that can be used to generate sampling directly through the python
interface. This is particularly useful when simulating simple systems, where
the overhead of starting up an external program is large compared to the actual
cost of computing the trajectory segment. Other use cases might include running
sampling with software that has a Python API (e.g. `OpenMM
<https://simtk.org/home/openmm>`_).

In order to create a custom propagator, users must define a class that inherits
from ``WESTPropagator`` and implement three methods:

- ``get_pcoord(self, state)``: Get the progress coordinate of the given basis
  or initial state.
- ``gen_istate(self, basis_state, initial_state)``: Generate a new initial
  state from the given basis state. This method is optional if ``gen_istates``
  is set to ``False`` in the propagation section of the configuration file,
  which is the default setting.
- ``propagate(self, segments)``: Propagate one or more segments, including any
  necessary per-iteration setup and teardown for this propagator.

There are also two stubs that that, if overridden, provide a mechanism for
modifying the simulation before or after the iteration:

- ``prepare_iteration(self, n_iter, segments)``: Perform any necessary
  per-iteration preparation. This is run by the work manager.
- ``finalize_iteration(self, n_iter, segments)``: Perform any necessary
  post-iteration cleanup. This is run by the work manager.

Several examples of custom propagators are available:

- `1D Over-damped Langevin dynamics
  <https://github.com/westpa/westpa/blob/master/lib/examples/odld/odld_system.py>`_
- `2D Langevin dynamics
  <https://bitbucket.org/joshua.adelman/stringmethodexamples/src/tip/examples/DicksonRingPotential/we_base/system.py>`_
- `Langevin dynamics - CA atom Elastic Network Model
  <https://bitbucket.org/joshua.adelman/stringmethodexamples/src/tip/examples/ElasticNetworkModel/we_base/system.py>`_

Configuration File
------------------

The configuration of a WESTPA simulation is specified using a plain text file
written in `YAML <http://en.wikipedia.org/wiki/YAML>`_. This file specifies,
among many other things, the length of the simulation, which modules should be
loaded for specifying the system, how external data should be organized on the
file system, and which plugins should used. YAML is a hierarchical format and
WESTPA organizes the configuration settings into blocks for each component.
While below, the configuration file will be referred to as **west.cfg**, the
user is free to name the configuration file something else. Most of the scripts
and tools that WESTPA provides, however, require that the name of the
configuration file be specified if the default name is not used.

The top most heading in *west.cfg* should be specified as:::

  ---
  west:
      ...

with all sub-section specified below it. A complete example can be found for
the NaCl example:
https://github.com/westpa/westpa/blob/master/lib/examples/nacl_gmx/west.cfg

In the following section, the specifications for each section of the file can
be found, along with default parameters and descriptions. Required parameters
are indicated as REQUIRED.::

  ---
  west:
      ...
      system:
          driver: REQUIRED 
          module_path: []  

The ``driver`` parameter must be set to a subclass of ``WESTSystem``, and given
in the form *module.class*. The ``module_path`` parameter is appended to the
system path and indicates where the class is defined.::

  ---
  west:
      ...
      we:
          adjust_counts: True
          weight_split_threshold: 2.0
          weight_merge_cutoff: 1.0

The ``we`` section section specifies parameters related to the Huber and Kim
resampling algorithm. WESTPA implements a variation of the method, in which
setting ``adust_counts`` to ``True`` strictly enforces that the number of
replicas per bin is exactly ``system.bin_target_counts``. Otherwise, the number
of replicas per is allowed to fluctuate as in the original implementation of
the algorithm. Adjusting the counts can improve load balancing for parallel
simulations. Replicas with weights greater than ``weight_split_threshold``
times the ideal weight per bin are tagged as candidates for splitting. Replicas
with weights less than ``weight_merge_cutoff`` times the ideal weight per bin
are candidates for merging.::

  ---
  west:
      ...
      propagation:
          gen_istates: False
          block_size: 1
          save_transition_matrices: False
          max_run_wallclock: None
          max_total_iterations: None

- ``gen_istates``: Boolean specifying whether to generate initial states from
  the basis states. The executable propagator defines a specific configuration
  block (*add internal link to other section*), and custom propagators should
  override the ``WESTPropagator.gen_istate()`` method.
- ``block_size``: An integer defining how many segments should be passed to a
  worker at a time. When using the serial work manager, this value should be
  set to the maximum number of segments per iteration to avoid significant
  overhead incurred by the locking mechanism in the WMFutures framework.
  Parallel work managers might benefit from setting this value greater than one
  in some instances to decrease network communication load.
- ``save_transition_matrices``:
- ``max_run_wallclock``: A time in dd:hh:mm:ss or hh:mm:ss specifying the
  maximum wallclock time of a particular WESTPA run. If running on a batch
  queuing system, this time should be set to less than the job allocation time
  to ensure that WESTPA shuts down cleanly.
- ``max_total_iterations``: An integer value specifying the number of
  iterations to run. This parameter is checked against the last completed
  iteration stored in the HDF5 file, not the number of iterations completed for
  a specific run. The default value of ``None`` only stops upon external
  termination of the code.::

    ---
    west:
        ...
        data:
            west_data_file: REQUIRED
            aux_compression_threshold: 1048576
            iter_prec: 8
            datasets:
                -name: REQUIRED
                 h5path: 
                 store: True
                 load: False
                 dtype: 
                 scaleoffset: None 
                 compression: None
                 chunks: None
            data_refs:
                segment: 
                basis_state:
                initial_state:

- ``west_data_file``: The name of the main HDF5 data storage file for the
  WESTPA simulation.
- ``aux_compression_threshold``: The threshold in bytes for compressing the
  auxiliary data in a dataset on an iteration-by-iteration basis.
- ``iter_prec``: The length of the iteration index with zero-padding. For the
  default value, iteration 1 would be specified as iter_00000001.
- ``datasets``:
- ``data_refs``:
- plugins
- executable

Environmental Variables
-----------------------

There are a number of environmental variables that can be set by the user in
order to configure a WESTPA simulation:

- WEST_ROOT: path to the base directory containing the WESTPA install
- WEST_SIM_ROOT: path to the base directory of the WESTPA simulation
- WEST_PYTHON: path to python executable to run the WESTPA simulation
- WEST_PYTHONPATH: path to any additional modules that WESTPA will require to
  run the simulation
- WEST_KERNPROF: path to ``kernprof.py`` script to perform line-by-line
  profiling of a WESTPA simulation (see `python line_profiler
  <http://pythonhosted.org/line_profiler>`__). This is only required for users
  who need to profile specific methods in a running WESTPA simulation.

Work manager related environmental variables:

- WM_WORK_MANAGER
- WM_N_WORKERS

WESTPA makes available to any script executed by it (e.g. **runseg.sh**), a
number of environmental variables that are set dynamically by the executable
propagator from the running simulation.

Programs executed for an iteration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following environment variables are passed to programs executed on a
per-iteration basis, notably pre-iteration and post-iteration scripts.

=================== =================== =======================================
Variable            Possible values     Function
=================== =================== =======================================
WEST_CURRENT_ITER   Integer >=1         Current iteration number
=================== =================== =======================================

Programs executed for a segment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following environment variables are passed to programs executed on a
per-segment basis, notably dynamics propagation.

=============================== =========================== ===================
Variable                        Possible values             Function
=============================== =========================== ===================
WEST_CURRENT_ITER               Integer >=1                 Current iteration
                                                            number
WEST_CURRENT_SEG_ID             Integer >=0                 Current segment ID
WEST_CURRENT_SEG_DATA_REF       String                      General-purpose
                                                            reference, based on
                                                            current segment
                                                            information,
                                                            configured in
                                                            west.cfg. Usually
                                                            used for storage
                                                            paths
WEST_CURRENT_SEG_INITPOINT_TYPE Enumeration:                Whether this
                                SEG_INITPOINT_CONTINUES,    segment continues a
                                SEG_INITPOINT_NEWTRAJ       previous trajectory
                                                            or initiates a new
                                                            one.
WEST_PARENT_ID                  Integer                     Segment ID of
                                                            parent segment.
                                                            Negative for
                                                            initial points.
WEST_PARENT_DATA_REF            String                      General purpose
                                                            reference, based on
                                                            parent segment
                                                            information,
                                                            configured in
                                                            west.cfg. Usually
                                                            used for storage
                                                            paths
WEST_PCOORD_RETURN              Filename                    Where progress
                                                            coordinate data
                                                            must be stored
WEST_RAND16                     Integer                     16-bit random
                                                            integer
WEST_RAND32                     Integer                     32-bit random
                                                            integer
WEST_RAND64                     Integer                     64-bit random
                                                            integer
WEST_RAND128                    Integer                     128-bit random
                                                            integer
WEST_RANDFLOAT                  Floating-point              Random number in
                                                            [0,1).
=============================== =========================== ===================

Additionally for any additional datasets specified in the configuration file,
WESTPA automatically provides ``WEST_X_RETURN``, where ``X`` is the uppercase
name of the dataset. For example if the configuration file contains the
following::

  data:
      ...
      datasets: # dataset storage options
        - name: energy

WESTPA would make ``WEST_ENERGY_RETURN`` available.

Programs executed for a single point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Programs used for creating initial states from basis states (``gen_istate.sh``)
or extracting progress coordinates from structures (e.g. ``get_pcoord.sh``) are
provided the following environment variables:

======================= =============== =============== =======================
Variable                Available for   Possible values Function
======================= =============== =============== =======================
WEST_STRUCT_DATA_REF    All             String          General-purpose
                        single-point                    reference, usually a
                        calculations                    pathname, associated
                                                        with the basis/initial
                                                        state.
WEST_BSTATE_ID          get_pcoord for  Integer >= 0    Basis state ID
                        basis state,                   
                        gen_istate
WEST_BSTATE_DATA_REF    get_pcoord for  String          Basis state data
                        basis state,                    reference
                        gen_istate
WEST_ISTATE_ID          get_pcoord for  Integer >= 0    Inital state ID
                        initial state,
                        gen_istate
WEST_ISTATE_DATA_REF    get_pcoord for  String          Initial state data
                        initial state,                  references, usually a
                        gen_istate                      pathname
WEST_PCOORD_RETURN      get_pcoord for  Pathname        Where progress
                        basis or                        coordinate data is 
                        initial state                   expected to be found
                                                        after execution
======================= =============== =============== =======================

Plugins
-------

WESTPA has a extensible plugin architecture that allows the user to manipulate
the simulation at specified points during an iteration.

-  Activating plugins in the config file
-  Plugin execution order/priority

Weighted Ensemble Algorithm (Resampling)
----------------------------------------
