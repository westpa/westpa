.. _w_pdist:

w_pdist
=======

usage::

 w_pdist [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
               [--max-queue-length MAX_QUEUE_LENGTH] [-W WEST_H5FILE] [--first-iter N_ITER]
               [--last-iter N_ITER] [-b BINEXPR] [-o OUTPUT] [-C] [--loose]
               [--construct-dataset CONSTRUCT_DATASET | --dsspecs DSSPEC [DSSPEC ...]]
               [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
               [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
               [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
               [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
               [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
               [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
               [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

Calculate time-resolved, multi-dimensional probability distributions of WE
datasets.

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
Histogram binning
-----------------------------------------------------------------------------

By default, histograms are constructed with 100 bins in each dimension. This
can be overridden by specifying -b/--bins, which accepts a number of different
kinds of arguments::

  a single integer N
    N uniformly spaced bins will be used in each dimension.

  a sequence of integers N1,N2,... (comma-separated)
    N1 uniformly spaced bins will be used for the first dimension, N2 for the
    second, and so on.

  a list of lists [[B11, B12, B13, ...], [B21, B22, B23, ...], ...]
    The bin boundaries B11, B12, B13, ... will be used for the first dimension,
    B21, B22, B23, ... for the second dimension, and so on. These bin
    boundaries need not be uniformly spaced. These expressions will be
    evaluated with Python's ``eval`` construct, with ``np`` available for
    use [e.g. to specify bins using np.arange()].

The first two forms (integer, list of integers) will trigger a scan of all
data in each dimension in order to determine the minimum and maximum values,
which may be very expensive for large datasets. This can be avoided by
explicitly providing bin boundaries using the list-of-lists form.

Note that these bins are *NOT* at all related to the bins used to drive WE
sampling.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file produced (specified by -o/--output, defaulting to "pdist.h5")
may be fed to plothist to generate plots (or appropriately processed text or
HDF5 files) from this data. In short, the following datasets are created::

  ``histograms``
    Normalized histograms. The first axis corresponds to iteration, and
    remaining axes correspond to dimensions of the input dataset.

  ``/binbounds_0``
    Vector of bin boundaries for the first (index 0) dimension. Additional
    datasets similarly named (/binbounds_1, /binbounds_2, ...) are created
    for additional dimensions.

  ``/midpoints_0``
    Vector of bin midpoints for the first (index 0) dimension. Additional
    datasets similarly named are created for additional dimensions.

  ``n_iter``
    Vector of iteration numbers corresponding to the stored histograms (i.e.
    the first axis of the ``histograms`` dataset).

-----------------------------------------------------------------------------
Subsequent processing
-----------------------------------------------------------------------------

The output generated by this program (-o/--output, default "pdist.h5") may be
plotted by the ``plothist`` program. See ``plothist --help`` for more
information.

-----------------------------------------------------------------------------
Parallelization
-----------------------------------------------------------------------------

This tool supports parallelized binning, including reading of input data.
Parallel processing is the default. For simple cases (reading pre-computed
input data, modest numbers of segments), serial processing (--serial) may be
more efficient.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------

optional arguments::

  -h, --help            show this help message and exit
  -b BINEXPR, --bins BINEXPR
                        Use BINEXPR for bins. This may be an integer, which will be used for each
                        dimension of the progress coordinate; a list of integers (formatted as
                        [n1,n2,...]) which will use n1 bins for the first dimension, n2 for the second
                        dimension, and so on; or a list of lists of boundaries (formatted as [[a1, a2,
                        ...], [b1, b2, ...], ... ]), which will use [a1, a2, ...] as bin boundaries for
                        the first dimension, [b1, b2, ...] as bin boundaries for the second dimension,
                        and so on. (Default: 100 bins in each dimension.)
  -o OUTPUT, --output OUTPUT
                        Store results in OUTPUT (default: pdist.h5).
  -C, --compress        Compress histograms. May make storage of higher-dimensional histograms more
                        tractable, at the (possible extreme) expense of increased analysis time.
                        (Default: no compression.)
  --loose               Ignore values that do not fall within bins. (Risky, as this can make buggy bin
                        boundaries appear as reasonable data. Only use if you are sure of your bin
                        boundary specification.)

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

WEST input data options:
  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

iteration range::

  --first-iter N_ITER   Begin analysis at iteration N_ITER (default: 1).
  --last-iter N_ITER    Conclude analysis with N_ITER, inclusive (default: last completed iteration).

input dataset options::

  --construct-dataset CONSTRUCT_DATASET
                        Use the given function (as in module.function) to extract source data. This
                        function will be called once per iteration as function(n_iter, iter_group) to
                        construct data for one iteration. Data returned must be indexable as
                        [seg_id][timepoint][dimension]
  --dsspecs DSSPEC [DSSPEC ...]
                        Construct probability distribution from one or more DSSPECs.

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