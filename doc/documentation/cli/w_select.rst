w_select
========

usage::

 w_select [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                [--max-queue-length MAX_QUEUE_LENGTH] [-W WEST_H5FILE] [--first-iter N_ITER]
                [--last-iter N_ITER] [-p MODULE.FUNCTION] [-v] [-a] [-o OUTPUT]
                [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
                [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
                [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
                [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
                [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
                [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
                [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

Select dynamics segments matching various criteria. This requires a
user-provided prediate function. By default, only matching segments are
stored. If the -a/--include-ancestors option is given, then matching segments
and their ancestors will be stored.

Predicate function
------------------

Segments are selected based on a predicate function, which must be callable
as ``predicate(n_iter, iter_group)`` and return a collection of segment IDs
matching the predicate in that iteration.

The predicate may be inverted by specifying the -v/--invert command-line
argument.

Output format
-------------

The output file (-o/--output, by default "select.h5") contains the following
datasets::

  ``/n_iter`` [iteration]
    *(Integer)* Iteration numbers for each entry in other datasets.

  ``/n_segs`` [iteration]
    *(Integer)* Number of segment IDs matching the predicate (or inverted
    predicate, if -v/--invert is specified) in the given iteration.

  ``/seg_ids`` [iteration][segment]
    *(Integer)* Matching segments in each iteration. For an iteration
    ``n_iter``, only the first ``n_iter`` entries are valid. For example,
    the full list of matching seg_ids in the first stored iteration is
    ``seg_ids[0][:n_segs[0]]``.

  ``/weights`` [iteration][segment]
    *(Floating-point)* Weights for each matching segment in ``/seg_ids``.

Command-line arguments
----------------------

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
                        Maximum number of tasks that can be queued. Useful to limit RAM use for tasks that
                        have very large requests/response. Default: no limit.

WEST input data options::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

iteration range::

  --first-iter N_ITER   Begin analysis at iteration N_ITER (default: 1).
  --last-iter N_ITER    Conclude analysis with N_ITER, inclusive (default: last completed iteration).

selection options::

  -p MODULE.FUNCTION, --predicate-function MODULE.FUNCTION
                        Use the given predicate function to match segments. This function should take an
                        iteration number and the HDF5 group corresponding to that iteration and return a
                        sequence of seg_ids matching the predicate, as in ``match_predicate(n_iter,
                        iter_group)``.
  -v, --invert          Invert the match predicate.
  -a, --include-ancestors
                        Include ancestors of matched segments in output.

output options:
  -o OUTPUT, --output OUTPUT
                        Write output to OUTPUT (default: select.h5).

parallelization options::

  --serial              run in serial mode
  --parallel            run in parallel mode (using processes)
  --work-manager WORK_MANAGER
                        use the given work manager for parallel task distribution. Available work managers
                        are ('serial', 'threads', 'processes', 'zmq'); default is 'serial'
  --n-workers N_WORKERS
                        Use up to N_WORKERS on this host, for work managers which support this option. Use
                        0 for a dedicated server. (Ignored by work managers which do not support this
                        option.)

options for ZeroMQ ("zmq") work manager (master or node)::

  --zmq-mode MODE       Operate as a master (server) or a node (workers/client). "server" is a deprecated
                        synonym for "master" and "client" is a deprecated synonym for "node".
  --zmq-comm-mode COMM_MODE
                        Use the given communication mode -- TCP or IPC (Unix-domain) -- sockets for
                        communication within a node. IPC (the default) may be more efficient but is not
                        available on (exceptionally rare) systems without node-local storage (e.g. /tmp);
                        on such systems, TCP may be used instead.
  --zmq-write-host-info INFO_FILE
                        Store hostname and port information needed to connect to this instance in
                        INFO_FILE. This allows the master and nodes assisting in coordinating the
                        communication of other nodes to choose ports randomly. Downstream nodes read this
                        file with --zmq-read-host-info and know where how to connect.
  --zmq-read-host-info INFO_FILE
                        Read hostname and port information needed to connect to the master (or other
                        coordinating node) from INFO_FILE. This allows the master and nodes assisting in
                        coordinating the communication of other nodes to choose ports randomly, writing
                        that information with --zmq-write-host-info for this instance to read.
  --zmq-upstream-rr-endpoint ENDPOINT
                        ZeroMQ endpoint to which to send request/response (task and result) traffic toward
                        the master.
  --zmq-upstream-ann-endpoint ENDPOINT
                        ZeroMQ endpoint on which to receive announcement (heartbeat and shutdown
                        notification) traffic from the master.
  --zmq-downstream-rr-endpoint ENDPOINT
                        ZeroMQ endpoint on which to listen for request/response (task and result) traffic
                        from subsidiary workers.
  --zmq-downstream-ann-endpoint ENDPOINT
                        ZeroMQ endpoint on which to send announcement (heartbeat and shutdown
                        notification) traffic toward workers.
  --zmq-master-heartbeat MASTER_HEARTBEAT
                        Every MASTER_HEARTBEAT seconds, the master announces its presence to workers.
  --zmq-worker-heartbeat WORKER_HEARTBEAT
                        Every WORKER_HEARTBEAT seconds, workers announce their presence to the master.
  --zmq-timeout-factor FACTOR
                        Scaling factor for heartbeat timeouts. If the master doesn't hear from a worker in
                        WORKER_HEARTBEAT*FACTOR, the worker is assumed to have crashed. If a worker
                        doesn't hear from the master in MASTER_HEARTBEAT*FACTOR seconds, the master is
                        assumed to have crashed. Both cases result in shutdown.
  --zmq-startup-timeout STARTUP_TIMEOUT
                        Amount of time (in seconds) to wait for communication between the master and at
                        least one worker. This may need to be changed on very large, heavily-loaded
                        computer systems that start all processes simultaneously.
  --zmq-shutdown-timeout SHUTDOWN_TIMEOUT
                        Amount of time (in seconds) to wait for workers to shut down.


westpa.cli.tools.w\_select module
---------------------------------

.. automodule:: westpa.cli.tools.w_select
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
