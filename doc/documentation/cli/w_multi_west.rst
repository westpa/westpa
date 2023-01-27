w_multi_west
============

The ``w_multi_west`` tool combines multiple WESTPA simulations into a single aggregate simulation to facilitate the analysis of the set of simulations. In particular, the tool creates a single ``west.h5`` file that contains all of the data from the west.h5 files of the individual simulations. Each iteration x in the new file contains all of the segments from iteration x from each of the set of simulation, all normalized to the total weight.

Overview
--------

usage::

 w_multi_west [-h] [-m master] [-n sims] [--quiet | --verbose | --debug] [--version] 
              [-W WEST_H5FILE] [-a aux] [--auxall] [--ibstates]
              [--serial | --parallel | --work-manager WORK_MANAGER] [--n-workers N_WORKERS]
              [--zmq-mode MODE] [--zmq-comm-mode COMM_MODE] [--zmq-write-host-info INFO_FILE]
              [--zmq-read-host-info INFO_FILE] [--zmq-upstream-rr-endpoint ENDPOINT]
              [--zmq-upstream-ann-endpoint ENDPOINT] [--zmq-downstream-rr-endpoint ENDPOINT]
              [--zmq-downstream-ann-endpoint ENDPOINT] [--zmq-master-heartbeat MASTER_HEARTBEAT]
              [--zmq-worker-heartbeat WORKER_HEARTBEAT] [--zmq-timeout-factor FACTOR]
              [--zmq-startup-timeout STARTUP_TIMEOUT] [--zmq-shutdown-timeout SHUTDOWN_TIMEOUT]

optional arguments::

   -h, --help           show this help message and exit

General options::
  -m, --master directory
                        Master path of simulations where all the smaller simulations are stored 
                        (default: Current Directory)
  -n, --sims n          Number of simulation directories. Assumes leading zeros. (default: 0)
  --quiet               emit only essential information
  --verbose             emit extra information
  --version             show program's version number and exit

 
Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for more
information on the general options.


Input/output options
~~~~~~~~~~~~~~~~~~~~

These arguments allow the user to specify where to read input simulation result
data and where to output calculated progress coordinate probability
distribution data.

Both input and output files are *hdf5* format::

  -W, --west, --WEST_H5FILE file
    The name of the main .h5 file inside each simulation directory. (Default: west.h5)

  -o, --output file
    Store this tool's output in file. (Default: multi.h5)

  -a, --aux auxdata
    Name of additional auxiliary dataset to be combined. Can be called multiple times. 
    (Default: None)

  -aa, --auxall
    Combine all auxiliary datsets as labeled in ``west.h5`` in folder 01. (Default: False)

  -nr, --no-reweight
    Do not perform reweighting. (Default: False)

  -ib, --ibstates
    Attempt to combine ``ibstates`` dataset if the basis states are identical across 
    all simulations. Needed when tracing with ``westpa.analysis``. (Default: False)

Examples
--------
If you have five simulations, set up your directory such that you 
have five directories are named numerically with leading zeroes, and each 
directory contains a ``west.h5`` file. For this example, each ``west.h5`` 
also contains an auxiliary dataset called ``RMSD``. If you run ``ls``, you 
will see the following output::

  01 02 03 04 05 

To run the w_multi_west tool, do the following::

  w_multi_west.py -m . -n 5 --aux=RMSD

If you used any custom WESTSystem, include that in the directory where you 
run the code.

To proceed in analyzing the aggregated simulation data as a single 
simulation, rename the output file ``multi.h5`` to ``west.h5``.

westpa.cli.tools.w\_multi\_west module
--------------------------------------

.. automodule:: westpa.cli.tools.w_multi_west
   :members:
   :undoc-members:
   :show-inheritance:
   :imported-members:
