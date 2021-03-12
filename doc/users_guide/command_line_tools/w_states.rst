.. _w_states:

w_states
========

usage::

 w_states [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
                [--show | --append | --replace] [--bstate-file BSTATE_FILE] [--bstate BSTATES]
                [--tstate-file TSTATE_FILE] [--tstate TSTATES]
                [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
                [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
                [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
                [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
                [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
                [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
                [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

Display or manipulate basis (initial) or target (recycling) states for a WEST simulation. By default,
states are displayed (or dumped to files). If ``--replace`` is specified, all basis/target states are
replaced for the next iteration. If ``--append`` is specified, the given target state(s) are appended to
the list for the next iteration. Appending basis states is not permitted, as this would require
renormalizing basis state probabilities in ways that may be error-prone. Instead, use ``w_states --show
--bstate-file=bstates.txt`` and then edit the resulting ``bstates.txt`` file to include the new desired
basis states, then use ``w_states --replace --bstate-file=bstates.txt`` to update the WEST HDF5 file
appropriately.

optional arguments::

  -h, --help            show this help message and exit
  --bstate-file BSTATE_FILE
                        Read (--append/--replace) or write (--show) basis state names, probabilities, and
                        data references from/to BSTATE_FILE.
  --bstate BSTATES      Add the given basis state (specified as a string 'label,probability[,auxref]') to
                        the list of basis states (after those specified in --bstate-file, if any). This
                        argument may be specified more than once, in which case the given states are
                        appended in the order they are given on the command line.
  --tstate-file TSTATE_FILE
                        Read (--append/--replace) or write (--show) target state names and representative
                        progress coordinates from/to TSTATE_FILE
  --tstate TSTATES      Add the given target state (specified as a string 'label,pcoord0[,pcoord1[,...]]')
                        to the list of target states (after those specified in the file given by
                        --tstates-from, if any). This argument may be specified more than once, in which
                        case the given states are appended in the order they appear on the command line.

general options::

  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --debug               enable extra checks and emit copious information
  --version             show program's version number and exit

modes of operation::

  --show                Display current basis/target states (or dump to files).
  --append              Append the given basis/target states to those currently in use.
  --replace             Replace current basis/target states with those specified.

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