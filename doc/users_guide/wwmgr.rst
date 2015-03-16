.. _wwmgr:

WEST Work Manager
=================

Introduction
------------

WWMGR is the parallel task distribution framework originally included as part
of the WEMD source. It was extracted to permit independent development, and
(more importantly) independent testing. A number of different schemes can be
selected at run-time for distributing work across multiple cores/nodes, as
follows:

=========== =================================================== =========== =========== ===============================
Name        Implementation                                      Multi-Core  Multi-Node  Appropriate For
=========== =================================================== =========== =========== ===============================
serial      None                                                No          No          Testing, minimizing overhead
                                                                                        when dynamics is inexpensive
threads     Python "threading" module                           Yes         No          Dynamics propagated by external
                                                                                        executables, large amounts of
                                                                                        data transferred per segment
processes   Python "multiprocessing" module                     Yes         No          Dynamics propagated by Python
                                                                                        routines, modest amounts of
                                                                                        data transferred per segment
mpi         `mpi4py <http://mpi4py.scipy.org/>`_                Yes         Yes         Distributing calculations
            compiled and linked against system MPI                                      across multiple nodes. Start
                                                                                        with this on your cluster of
                                                                                        choice. 
zmq         `ZeroMQ <http://www.zeromq.org/>`_                  Yes         Yes         Distributing calculations
            and `PyZMQ <http://zeromq.github.com/pyzmq/>`_                              across multiple nodes. Use this
                                                                                        if MPI does not work properly
                                                                                        on your cluster (particularly
                                                                                        for spawning child processes).
=========== =================================================== =========== =========== ===============================

Environment variables
---------------------

For controlling task distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the original WEMD work managers were controlled by command-line options
and entries in wemd.cfg, the new work manager is controlled using command-line
options or environment variables (much like OpenMP). These variables are as
follow:

=========================== ======================= =================================== ===============================
Variable                    Applicable to           Default                             Meaning
=========================== ======================= =================================== ===============================
WM_WORK_MANAGER             (none)                  processes                           Use the given task distribution
                                                                                        system: "serial", "threads",
                                                                                        "processes", or "zmq"
WM_N_WORKERS                threads, processes, zmq number of cores in machine          Use this number of workers. In
                                                                                        the case of zmq, use this many
                                                                                        workers on the current machine
                                                                                        only (can be set independently
                                                                                        on different nodes).
WM_ZMQ_MODE                 zmq                     server                              Start as a server ("server") or
                                                                                        a client ("client"). Servers
                                                                                        coordinate a given calculation,
                                                                                        and clients execute tasks
                                                                                        related to that calculation.
WM_ZMQ_TASK_TIMEOUT         zmq                     60                                  Time (in seconds) after which a
                                                                                        worker will be considered hung,
                                                                                        terminated, and restarted. This
                                                                                        **must** be updated for
                                                                                        long-running dynamics segments.
                                                                                        Set to zero to disable hang
                                                                                        checks entirely.
WM_ZMQ_TASK_ENDPOINT        zmq                     Random port                         Master distributes tasks at
                                                                                        this address
WM_ZMQ_RESULT_ENDPOINT      zmq                     Random port                         Master receives task results at
                                                                                        this address                                                                                                                                                           |
WM_ZMQ_ANNOUNCE_ENDPOINT    zmq                     Random port                         Master publishes announcements
                                                                                        (such as "shut down now") at
                                                                                        this address
WM_ZMQ_SERVER_INFO          zmq                     ``zmq_server_info_PID_ID.json``     A file describing the above
                                                    (where PID is a process ID and      endpoints can be found here (to
                                                    ID is a nearly random hex number)   ease cluster-wide startup)
=========================== ======================= =================================== ===============================

For passing information to workers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One environment variable is made available by multi-process work managers
(processes and ZMQ) to help clients configure themselves (e.g. select an
appropriate GPU on a multi-GPU node):

=============== =============== ===============================================
Variable        Applicable to   Meaning
=============== =============== ===============================================
WM_PROCESS_ID   processes, zmq  Contains an integer, 0 based, identifying the
                                process among the set of processes started on a
                                given node.
=============== =============== ===============================================

The ZeroMQ work manager for clusters
------------------------------------

The ZeroMQ ("zmq") work manager can be used for both single-machine and
cluster-wide communication. Communication occurs over sockets using the `ZeroMQ
<http://www.zeromq.org/>`_ messaging protocol. Within nodes, `Unix sockets
<http://en.wikipedia.org/wiki/UNIX_socket>`_ are used for efficient
communication, while between nodes, TCP sockets are used. This also minimizes
the number of open sockets on the master node.

The quick and dirty guide to using this on a cluster is as follows::

    source env.sh
    export WM_WORK_MANAGER=zmq
    export WM_ZMQ_COMM_MODE=tcp
    export WM_ZMQ_SERVER_INFO=$WEST_SIM_ROOT/wemd_server_info.json

    w_run & 

    # manually run w_run on each client node, as appropriate for your batch system
    # e.g. qrsh -inherit for Grid Engine, or maybe just simple SSH

    for host in $(cat $TMPDIR/machines | sort | uniq); do
       qrsh -inherit -V $host $PWD/node-ltc1.sh &
    done
