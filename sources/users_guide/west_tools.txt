WEST Tools
==========

The command line tools included with the WESTPA software package are broadly
separable into two categories: **Tools for initializing a simulation** and
**tools for analyzing results**.

Command function can be user defined and modified. The particular parameters of
different command line tools are specified, in order of precedence, by:

- User specified command line arguments
- User defined environmental variables
- Package defaults

This page focuses on outlining the general functionality of the command line
tools and providing an overview of command line arguments that are shared by
multiple tools. See the :ref:`index of command-line tools
<command_line_tool_index>` for a more comprehensive overview of each tool.

Overview
--------

All tools are located in the ``$WEST_ROOT/bin`` directory, where the shell
variable ``WEST_ROOT`` points to the path where the WESTPA package is located
on your machine.

You may wish to set this variable automatically by adding the following to your
``~/.bashrc`` or ``~/.profile`` file::

  export WEST_ROOT="$HOME/westpa"

where the path to the westpa suite is modified accordingly.

Tools for setting up and running a simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following commands to initialize, configure, and run a weighted
ensemble simulation. Command line arguments or environmental variables can be
set to specify the work managers for running the simulation, where
configuration data is read from, and the *HDF5* file in which results are
stored.

=================== ===========================================================
Command             Function 
=================== ===========================================================
:ref:`w_init`       Initializes simulation configuration files and environment.
                    Always run this command before starting a new simulation.
:ref:`w_bins`       Set up binning, progress coordinate 
:ref:`w_run`        Launches a simulation. Command arguments/environmental
                    variables can be included to specify the work managers and
                    simulation parameters
:ref:`w_truncate`   Truncates the weighted ensemble simulation from a given
                    iteration. 
=================== ===========================================================

Tools for analyzing simulation results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following command line tools are provided for analysis after running a
weighted ensemble simulation (and collecting the results in an HDF5 file).

With the exception of the plotting tool ``plothist``, all analysis tools read
from and write to *HDF5* type files.

=================== ===========================================================
Command             Function
=================== ===========================================================
:ref:`w_assign`     Assign walkers to bins and macrostates (using simulation
                    output as input). Must be done before some other analysis
                    tools (e.g. :ref:`w_kinetics`, :ref:`w_kinavg`)
:ref:`w_trace`      Trace the path of a given walker segment over a
                    user-specified number of simulation iterations.
:ref:`w_fluxanl`    Calculate average probability flux into user-defined
                    'target' state with relevant statistics. 
:ref:`w_pdist`      Construct a probability distribution of results (e.g.
                    progress coordinate membership) for subsequent plotting
                    with :ref:`plothist`.
:ref:`plothist`     Tool to plot output from other analysis tools (e.g.
                    :ref:`w_pdist`). 
=================== ===========================================================

General Command Line Options
----------------------------

The following arguments are shared by all command line tools::

  -r config file, --rcfile config file
    Use config file as the configuration file (Default: File named west.cfg)
  --quiet, --verbose, --debug
    Specify command tool output verbosity (Default: 'quiet' mode)
  --version
    Print WESTPA version number and exit
  -h, --help
    Output the help information for this command line tool and exit

A note on specifying a configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A *configuration file*, which should be stored in your simulation root
directory, is read by all command line tools. The *configuration file*
specifies parameters for general simulation setup, as well as the *hdf5* file
name where simulation data is stored and read by analysis tools.

If not specified, the **default configuration file** is assumed to be named
**west.cfg**.

You can override this to use configuration file *file* by either:

-  Setting the environmental variable ``WESTRC`` equal to *file*::

    export WESTRC=/path/to/westrcfile

-  Including the command line argument ``-r /path/to/westrcfile``

Work Manager Options
--------------------

Note: See :ref:`wwmgr overview <wwmgr>` for a more detailed explanation of the
work manager framework.

Work managers a used by a number of command-line tools to process more complex
tasks, especially in setting up and running simulations (i.e. :ref:`w_init` and
:ref:`w_run`) - in general, work managers are involved in tasks that require
multiprocessing and/or tasks distributed over multiple nodes in a cluster.

Overview
~~~~~~~~

The following command-line tools make use of work managers:

- :ref:`w_init`
- :ref:`w_run`

General work manager options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following are general options used for specifying the type of work
manager and number of cores::

  --wm-work-manager work_manager
    Specify which type of work manager to use, where the possible choices for
    work_manager are: {processes, gcserial, threads, mpi, or zmq}. See the
    wwmgr overview page <wwmgr>_ for more information on the different types of
    work managers (Default: gcprocesses)
  --wm-n-workers n_workers
    Specify the number of cores to use as gcn_workers, if the work manager you
    selected supports this option (work managers that do not will ignore this
    option). If using an gcmpi or zmq work manager, specify gc--wm-n-workers=0
    for a dedicated server (Default: Number of cores available on machine)

The ``mpi`` work manager is generally sufficient for most tasks that make use
of multiple nodes on a cluster. The ``zmq`` work manager is preferable if the
``mpi`` work manager does not work properly on your cluster or if you prefer to
have more explicit control over the distribution of communication tasks on your
cluster.

ZeroMQ ('zmq') work manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ZeroMQ work manager offers a number of additional options (all of
which are optional and have default values). All of these options focus
on whether the zmq work manager is set up as a server (i.e. task
distributor/ventilator) or client (task processor)::

  --wm-zmq-mode mode
    Options: {server or client}. Specify whether the ZMQ work manager on this
    node will operate as a server or a client (Default: server)

  --wm-zmq-info-file info_file 
    Specify the name of a temporary file to write (as a server) or read (as a
    client) socket connection endpoints (Default: server_x.json, where x is a
    unique identifier string)

  --wm-zmq-task-endpoint task_endpoint 
    Explicitly use task_endpoint to bind to (as server) or connect to (as
    client) for task distribution (Default: A randomly determined endpoint that
    is written or read from the specified info_file)

  --wm-zmq-result-endpoint result_endpoint 
    Explicitly use result_endpoint to bind to (as server) or connect to (as
    client) to distribute and collect task results (Default: A randomly
    determined endpoint that is written to or read from the specified
    info_file)

  --wm-zmq-announce-endpoint announce_endpoint 
    Explicitly use announce_endpoint to bind to (as server) or connect to (as
    client) to distribute central announcements (Default: A randomly determined
    endpoint that is written to or read from the specified info_file)

  --wm-zmq-heartbeat-interval interval 
    If a server, send an Im alive ping to connected clients every interval
    seconds; If a client, expect to hear a server ping every approximately
    interval seconds, or else assume the server has crashed and shutdown
    (Default: 600 seconds)

  --wm-zmq-task-timeout timeout 
    Kill worker processes/jobs after that take longer than timeout seconds to
    complete (Default: no time limit)

  --wm-zmq-client-comm-mode mode 
    Use the communication mode, mode, (options: {ipc for Unix sockets, or tcp
    for TCP/IP sockets}) to communicate with worker processes (Default: ipc)

Initializing/Running Simulations
--------------------------------

For a more complete overview of all the files necessary for setting up a
simulation, see the :ref:`user guide for setting up a simulation <setup>`
