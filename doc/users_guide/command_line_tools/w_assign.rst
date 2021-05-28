.. _w_assign:

w_assign
========

usage::

 w_assign [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                [--max-queue-length MAX_QUEUE_LENGTH] [-W WEST_H5FILE]
                [--bins-from-system | --bins-from-expr BINS_FROM_EXPR | --bins-from-function BINS_FROM_FUNCTION | --bins-from-file BINFILE | --bins-from-h5file]
                [--construct-dataset CONSTRUCT_DATASET | --dsspecs DSSPEC [DSSPEC ...]]
                [--states STATEDEF [STATEDEF ...] | --states-from-file STATEFILE |
                --states-from-function STATEFUNC] [-o OUTPUT] [--subsample] [--config-from-file]
                [--scheme-name SCHEME] [--serial | --parallel | --work-manager WORK_MANAGER]
                [--n-workers N_WORKERS] [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE]
                [--zmq-write-host-info INFO_FILE] [--zmq-read-host-info INFO_FILE]
                [--zmq-upstream-rr-endpoint ENDPOINT] [--zmq-upstream-ann-endpoint ENDPOINT]
                [--zmq-downstream-rr-endpoint ENDPOINT] [--zmq-downstream-ann-endpoint ENDPOINT]
                [--zmq-master-heartbeat MASTER_HEARTBEAT] [--zmq-worker-heartbeat WORKER_HEARTBEAT]
                [--zmq-timeout-factor FACTOR] [--zmq-startup-timeout STARTUP_TIMEOUT]
                [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

Assign walkers to bins, producing a file (by default named "assign.h5")
which can be used in subsequent analysis.

For consistency in subsequent analysis operations, the entire dataset
must be assigned, even if only a subset of the data will be used. This
ensures that analyses that rely on tracing trajectories always know the
originating bin of each trajectory.

-----------------------------------------------------------------------------
Source data
-----------------------------------------------------------------------------

Source data is provided either by a user-specified function
(--construct-dataset) or a list of "data set specifications" (--dsspecs).
If neither is provided, the progress coordinate dataset ''pcoord'' is used.

To use a custom function to extract or calculate data whose probability
distribution will be calculated, specify the function in standard Python
MODULE.FUNCTION syntax as the argument to --construct-dataset. This function
will be called as function(n_iter,iter_group), where n_iter is the iteration
whose data are being considered and iter_group is the corresponding group
in the main WEST HDF5 file (west.h5). The function must return data which can
be indexed as [segment][timepoint][dimension].

To use a list of data set specifications, specify --dsspecs and then list the
desired datasets one-by-one (space-separated in most shells). These data set
specifications are formatted as NAME[,file=FILENAME,slice=SLICE], which will
use the dataset called NAME in the HDF5 file FILENAME (defaulting to the main
WEST HDF5 file west.h5), and slice it with the Python slice expression SLICE
(as in [0:2] to select the first two elements of the first axis of the
dataset). The ``slice`` option is most useful for selecting one column (or
more) from a multi-column dataset, such as arises when using a progress
coordinate of multiple dimensions.

-----------------------------------------------------------------------------
Specifying macrostates
-----------------------------------------------------------------------------

Optionally, kinetic macrostates may be defined in terms of sets of bins.
Each trajectory will be labeled with the kinetic macrostate it was most
recently in at each timepoint, for use in subsequent kinetic analysis.
This is required for all kinetics analysis (w_kintrace and w_kinmat).

There are three ways to specify macrostates:

  1. States corresponding to single bins may be identified on the command
     line using the --states option, which takes multiple arguments, one for
     each state (separated by spaces in most shells). Each state is specified
     as a coordinate tuple, with an optional label prepended, as in
     ``bound:1.0`` or ``unbound:(2.5,2.5)``. Unlabeled states are named
     ``stateN``, where N is the (zero-based) position in the list of states
     supplied to --states.

  2. States corresponding to multiple bins may use a YAML input file specified
     with --states-from-file. This file defines a list of states, each with a
     name and a list of coordinate tuples; bins containing these coordinates
     will be mapped to the containing state. For instance, the following
     file::

        ---
        states:
          - label: unbound
            coords:
              - [9.0, 1.0]
              - [9.0, 2.0]
          - label: bound
            coords:
              - [0.1, 0.0]

     produces two macrostates: the first state is called "unbound" and
     consists of bins containing the (2-dimensional) progress coordinate
     values (9.0, 1.0) and (9.0, 2.0); the second state is called "bound"
     and consists of the single bin containing the point (0.1, 0.0).

  3. Arbitrary state definitions may be supplied by a user-defined function,
     specified as --states-from-function=MODULE.FUNCTION. This function is
     called with the bin mapper as an argument (``function(mapper)``) and must
     return a list of dictionaries, one per state. Each dictionary must contain
     a vector of coordinate tuples with key "coords"; the bins into which each
     of these tuples falls define the state. An optional name for the state
     (with key "label") may also be provided.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "assign.h5") contains the following
attributes datasets::

  ``nbins`` attribute
    *(Integer)* Number of valid bins. Bin assignments range from 0 to
    *nbins*-1, inclusive.

  ``nstates`` attribute
    *(Integer)* Number of valid macrostates (may be zero if no such states are
    specified). Trajectory ensemble assignments range from 0 to *nstates*-1,
    inclusive, when states are defined.

  ``/assignments`` [iteration][segment][timepoint]
    *(Integer)* Per-segment and -timepoint assignments (bin indices).

  ``/npts`` [iteration]
    *(Integer)* Number of timepoints in each iteration.

  ``/nsegs`` [iteration]
    *(Integer)* Number of segments in each iteration.

  ``/labeled_populations`` [iterations][state][bin]
    *(Floating-point)* Per-iteration and -timepoint bin populations, labeled
    by most recently visited macrostate. The last state entry (*nstates-1*)
    corresponds to trajectories initiated outside of a defined macrostate.

  ``/bin_labels`` [bin]
    *(String)* Text labels of bins.

When macrostate assignments are given, the following additional datasets are
present::

  ``/trajlabels`` [iteration][segment][timepoint]
    *(Integer)* Per-segment and -timepoint trajectory labels, indicating the
    macrostate which each trajectory last visited.

  ``/state_labels`` [state]
    *(String)* Labels of states.

  ``/state_map`` [bin]
    *(Integer)* Mapping of bin index to the macrostate containing that bin.
    An entry will contain *nbins+1* if that bin does not fall into a
    macrostate.

Datasets indexed by state and bin contain one more entry than the number of
valid states or bins. For *N* bins, axes indexed by bin are of size *N+1*, and
entry *N* (0-based indexing) corresponds to a walker outside of the defined bin
space (which will cause most mappers to raise an error). More importantly, for
*M* states (including the case *M=0* where no states are specified), axes
indexed by state are of size *M+1* and entry *M* refers to trajectories
initiated in a region not corresponding to a defined macrostate.

Thus, ``labeled_populations[:,:,:].sum(axis=1)[:,:-1]`` gives overall per-bin
populations, for all defined bins and
``labeled_populations[:,:,:].sum(axis=2)[:,:-1]`` gives overall
per-trajectory-ensemble populations for all defined states.

-----------------------------------------------------------------------------
Parallelization
-----------------------------------------------------------------------------

This tool supports parallelized binning, including reading/calculating input
data.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------

optional arguments::

  -h, --help            show this help message and exit

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

parallelization options::

  --max-queue-length MAX_QUEUE_LENGTH
                        Maximum number of tasks that can be queued. Useful to limit RAM use for tasks
                        that have very large requests/response. Default: no limit.

WEST input data options::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

binning options:
  --bins-from-system    Bins are constructed by the system driver specified in the WEST configuration
                        file (default where stored bin definitions not available).
  --bins-from-expr BINS_FROM_EXPR, --binbounds BINS_FROM_EXPR
                        Construct bins on a rectilinear grid according to the given BINEXPR. This must
                        be a list of lists of bin boundaries (one list of bin boundaries for each
                        dimension of the progress coordinate), formatted as a Python expression. E.g.
                        "[[0,1,2,4,inf],[-inf,0,inf]]". The numpy module and the special symbol "inf"
                        (for floating-point infinity) are available for use within BINEXPR.
  --bins-from-function BINS_FROM_FUNCTION, --binfunc BINS_FROM_FUNCTION
                        Supply an external function which, when called, returns a properly constructed
                        bin mapper which will then be used for bin assignments. This should be
                        formatted as "[PATH:]MODULE.FUNC", where the function FUNC in module MODULE
                        will be used; the optional PATH will be prepended to the module search path
                        when loading MODULE.
  --bins-from-file BINFILE, --binfile BINFILE
                        Load bin specification from the YAML file BINFILE. This currently takes the
                        form {'bins': {'type': 'RectilinearBinMapper', 'boundaries': [[boundset1],
                        [boundset2], ... ]}}; only rectilinear bin bounds are supported.
  --bins-from-h5file    Load bin specification from the data file being examined (default where stored
                        bin definitions available).

input dataset options::

  --construct-dataset CONSTRUCT_DATASET
                        Use the given function (as in module.function) to extract source data. This
                        function will be called once per iteration as function(n_iter, iter_group) to
                        construct data for one iteration. Data returned must be indexable as
                        [seg_id][timepoint][dimension]
  --dsspecs DSSPEC [DSSPEC ...]
                        Construct source data from one or more DSSPECs.

macrostate definitions::

  --states STATEDEF [STATEDEF ...]
                        Single-bin kinetic macrostate, specified by a coordinate tuple (e.g. '1.0' or
                        '[1.0,1.0]'), optionally labeled (e.g. 'bound:[1.0,1.0]'). States corresponding
                        to multiple bins must be specified with --states-from-file.
  --states-from-file STATEFILE
                        Load kinetic macrostates from the YAML file STATEFILE. See description above
                        for the appropriate structure.
  --states-from-function STATEFUNC
                        Load kinetic macrostates from the function STATEFUNC, specified as
                        module_name.func_name. This function is called with the bin mapper as an
                        argument, and must return a list of dictionaries {'label': state_label,
                        'coords': 2d_array_like} one for each macrostate; the 'coords' entry must
                        contain enough rows to identify all bins in the macrostate.

other options::

  -o OUTPUT, --output OUTPUT
                        Store results in OUTPUT (default: assign.h5).
  --subsample           Determines whether or not the data should be subsampled. This is rather useful
                        for analysing steady state simulations.
  --config-from-file    Load bins/macrostates from a scheme specified in west.cfg.
  --scheme-name SCHEME  Name of scheme specified in west.cfg.

parallelization options::

  --serial              run in serial mode
  --parallel            run in parallel mode (using processes)
  --work-manager WORK_MANAGER
                        use the given work manager for parallel task distribution. Available work
                        managers are ('serial', 'threads', 'processes', 'zmq'); default is 'processes'
  --n-workers N_WORKERS
                        Use up to N_WORKERS on this host, for work managers which support this option.
                        Use 0 for a dedicated server. (Ignored by work managers which do not support
                        this option.)

options for ZeroMQ ("zmq") work manager (master or node):
  --zmq-mode MODE       Operate as a master (server) or a node (workers/client). "server" is a
                        deprecated synonym for "master" and "client" is a deprecated synonym for
                        "node".
  --zmq-comm-mode COMM_MODE
                        Use the given communication mode -- TCP or IPC (Unix-domain) -- sockets for
                        communication within a node. IPC (the default) may be more efficient but is not
                        available on (exceptionally rare) systems without node-local storage (e.g.
                        /tmp); on such systems, TCP may be used instead.
  --zmq-write-host-info INFO_FILE
                        Store hostname and port information needed to connect to this instance in
                        INFO_FILE. This allows the master and nodes assisting in coordinating the
                        communication of other nodes to choose ports randomly. Downstream nodes read
                        this file with --zmq-read-host-info and know where how to connect.
  --zmq-read-host-info INFO_FILE
                        Read hostname and port information needed to connect to the master (or other
                        coordinating node) from INFO_FILE. This allows the master and nodes assisting
                        in coordinating the communication of other nodes to choose ports randomly,
                        writing that information with --zmq-write-host-info for this instance to read.
  --zmq-upstream-rr-endpoint ENDPOINT
                        ZeroMQ endpoint to which to send request/response (task and result) traffic
                        toward the master.
  --zmq-upstream-ann-endpoint ENDPOINT
                        ZeroMQ endpoint on which to receive announcement (heartbeat and shutdown
                        notification) traffic from the master.
  --zmq-downstream-rr-endpoint ENDPOINT
                        ZeroMQ endpoint on which to listen for request/response (task and result)
                        traffic from subsidiary workers.
  --zmq-downstream-ann-endpoint ENDPOINT
                        ZeroMQ endpoint on which to send announcement (heartbeat and shutdown
                        notification) traffic toward workers.
  --zmq-master-heartbeat MASTER_HEARTBEAT
                        Every MASTER_HEARTBEAT seconds, the master announces its presence to workers.
  --zmq-worker-heartbeat WORKER_HEARTBEAT
                        Every WORKER_HEARTBEAT seconds, workers announce their presence to the master.
  --zmq-timeout-factor FACTOR
                        Scaling factor for heartbeat timeouts. If the master doesn't hear from a worker
                        in WORKER_HEARTBEAT*FACTOR, the worker is assumed to have crashed. If a worker
                        doesn't hear from the master in MASTER_HEARTBEAT*FACTOR seconds, the master is
                        assumed to have crashed. Both cases result in shutdown.
  --zmq-startup-timeout STARTUP_TIMEOUT
                        Amount of time (in seconds) to wait for communication between the master and at
                        least one worker. This may need to be changed on very large, heavily-loaded
                        computer systems that start all processes simultaneously.
  --zmq-shutdown-timeout SHUTDOWN_TIMEOUT
                        Amount of time (in seconds) to wait for workers to shut down.