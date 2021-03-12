.. _w_crawl:

w_crawl
=======

usage::

 w_crawl [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
               [--max-queue-length MAX_QUEUE_LENGTH] [-W WEST_H5FILE] [--first-iter N_ITER]
               [--last-iter N_ITER] [-c CRAWLER_INSTANCE]
               [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
               [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
               [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
               [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
               [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
               [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
               [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]
               task_callable

Crawl a weighted ensemble dataset, executing a function for each iteration.
This can be used for postprocessing of trajectories, cleanup of datasets,
or anything else that can be expressed as "do X for iteration N, then do
something with the result". Tasks are parallelized by iteration, and
no guarantees are made about evaluation order.

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

iteration range::

  --first-iter N_ITER   Begin analysis at iteration N_ITER (default: 1).
  --last-iter N_ITER    Conclude analysis with N_ITER, inclusive (default: last completed iteration).

task options::

  -c CRAWLER_INSTANCE, --crawler-instance CRAWLER_INSTANCE
                        Use CRAWLER_INSTANCE (specified as module.instance) as an instance of
                        WESTPACrawler to coordinate the calculation. Required only if initialization,
                        finalization, or task result processing is required.
  task_callable         Run TASK_CALLABLE (specified as module.function) on each iteration. Required.

parallelization options::

  --serial              run in serial mode
  --parallel            run in parallel mode (using processes)
  --work-manager WORK_MANAGER
                        use the given work manager for parallel task distribution. Available work
                        managers are ('serial', 'threads', 'processes', 'zmq'); default is 'serial'
  --n-workers N_WORKERS
                        Use up to N_WORKERS on this host, for work managers which support this option.
                        Use 0 for a dedicated server. (Ignored by work managers which do not support
                        this option.)

options for ZeroMQ ("zmq") work manager (master or node)::

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