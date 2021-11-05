w_red
=====

usage::

 w_red [-h] [-r RCFILE] [--quiet] [--verbose] [--version] [--max-queue-length MAX_QUEUE_LENGTH]
             [--debug] [--terminal]
             [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
             [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
             [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
             [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
             [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
             [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
             [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

optional arguments::

  -h, --help            show this help message and exit

general options:
  -r RCFILE, --rcfile RCFILE
                        use RCFILE as the WEST run-time configuration file (default: west.cfg)
  --quiet               emit only essential information
  --verbose             emit extra information
  --version             show program's version number and exit

parallelization options::

  --max-queue-length MAX_QUEUE_LENGTH
                        Maximum number of tasks that can be queued. Useful to limit RAM use for tasks that
                        have very large requests/response. Default: no limit.

parallelization options::

  --serial              run in serial mode
  --parallel            run in parallel mode (using processes)
  --work-manager WORK_MANAGER

  
westpa.cli.tools.w\_red module
------------------------------

.. automodule:: westpa.cli.tools.w_red
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
