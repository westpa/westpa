.. _w_assign:

w_assign
========

``w_assign`` uses simulation output to assign walkers to user-specified bins
and macrostates. These assignments are required for some other simulation
tools, namely ``w_kinetics`` and ``w_kinavg``.

``w_assign`` supports parallelization (see `general work manager options
<UserGuide:ToolRefs#Work_Manager_Options>`_ for more on command line options
to specify a work manager).

Overview
--------

Usage::

  w_assign [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                 [-W WEST_H5FILE] [-o OUTPUT]
                 [--bins-from-system | --bins-from-expr BINS_FROM_EXPR | --bins-from-function BINS_FROM_FUNCTION]
                 [-p MODULE.FUNCTION]
                 [--states STATEDEF [STATEDEF ...] | --states-from-file STATEFILE | --states-from-function STATEFUNC]
                 [--wm-work-manager WORK_MANAGER] [--wm-n-workers N_WORKERS]
                 [--wm-zmq-mode MODE] [--wm-zmq-info INFO_FILE]
                 [--wm-zmq-task-endpoint TASK_ENDPOINT]
                 [--wm-zmq-result-endpoint RESULT_ENDPOINT]
                 [--wm-zmq-announce-endpoint ANNOUNCE_ENDPOINT]
                 [--wm-zmq-listen-endpoint ANNOUNCE_ENDPOINT]
                 [--wm-zmq-heartbeat-interval INTERVAL]
                 [--wm-zmq-task-timeout TIMEOUT]
                 [--wm-zmq-client-comm-mode MODE]

Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for
more information on the general options.

Input/output Options
--------------------

::

  -W, --west-data /path/to/file

      Read simulation result data from file *file*. (**Default:** The
      *hdf5* file specified in the configuration file, by default
      **west.h5**)
  
  -o, --output /path/to/file
      Write assignment results to file *outfile*. (**Default:** *hdf5*
      file **assign.h5**)

Binning Options
---------------

Specify how binning is to be assigned to the dataset.::

  --bins-from-system
    Use binning scheme specified by the system driver; system driver can be
    found in the west configuration file, by default named **west.cfg**
    (**Default binning**)

  --bins-from-expr bin_expr
    Use binning scheme specified in *``bin_expr``*, which takes the form a
    Python list of lists, where each inner list corresponds to the binning a
    given dimension. (for example, "[[0,1,2,4,inf],[-inf,0,inf]]" specifies bin
    boundaries for two dimensional progress coordinate. Note that this option
    accepts the special symbol 'inf' for floating point infinity

  --bins-from-function bin_func
    Bins specified by calling an external function *``bin_func``*.
    *``bin_func``* should be formatted as '[PATH:]module.function', where the
    function 'function' in module 'module' will be used

Macrostate Options
------------------

You can optionally specify how to assign user-defined macrostates. Note
that macrostates must be assigned for subsequent analysis tools, namely
``w_kinetics`` and ``w_kinavg``.::

  --states statedef [statedef ...]
    Specify a macrostate for a single bin as *``statedef``*, formatted
    as a coordinate tuple where each coordinate specifies the bin to
    which it belongs, for instance:
    '[1.0, 2.0]' assigns a macrostate corresponding to the bin that
    contains the (two-dimensional) progress coordinates 1.0 and 2.0.
    Note that a macrostate label can optionally by specified, for
    instance: 'bound:[1.0, 2.0]' assigns the corresponding bin
    containing the given coordinates the macrostate named 'bound'. Note
    that multiple assignments can be specified with this command, but
    only one macrostate per bin is possible - if you wish to specify
    multiple bins in a single macrostate, use the
    *``--states-from-file``* option.

  --states-from-file statefile
    Read macrostate assignments from *yaml* file *``statefile``*. This
    option allows you to assign multiple bins to a single macrostate.
    The following example shows the contents of *``statefile``* that
    specify two macrostates, bound and unbound, over multiple bins with
    a two-dimensional progress coordinate:

  ---
  states:
    - label: unbound
      coords:
        - [9.0, 1.0]
        - [9.0, 2.0]
    - label: bound
      coords:
        - [0.1, 0.0]

Specifying Progress Coordinate
------------------------------

By default, progress coordinate information for each iteration is taken from
*pcoord* dataset in the specified input file (which, by default is *west.h5*).
Optionally, you can specify a function to construct the progress coordinate for
each iteration - this may be useful to consolidate data from several sources or
otherwise preprocess the progress coordinate data.::

  --construct-pcoord module.function, -p module.function
    Use the function *module.function* to construct the progress
    coordinate for each iteration. This will be called once per
    iteration as *function(n_iter, iter_group)* and should return an
    array indexable as [seg_id][timepoint][dimension]. The
    **default** function returns the 'pcoord' dataset for that iteration
    (i.e. the function executes return iter_group['pcoord'][...])

Examples
--------
