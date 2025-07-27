w_timings
=========

``w_timings`` prints timing information for a WESTPA simulation.

Overview
--------

Usage:

.. code-block:: shell

  w_timings [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
            [-W WEST_H5FILE] [--first-iter N_ITER] [--last-iter N_ITER] [-t TAU]

Optional arguments:

.. code-block:: shell

  -t TAU, --tau TAU     WE resampling interval (format: <value>_<unit>, where
                        <value> is a positive integer and <unit> is 'as', 'fs',
                        'ps', 'ns', 'us', 'ms', 's', 'm', 'h', 'D', or 'W').

Examples
--------

The default output includes the wall-clock time:

.. code-block:: console

  $ w_timings -W west.h5
  Iterations: 50
  Total segments: 9985
  Wall-clock time: 0:00:01.452625

To determine the simulated time, the resampling time (``-t`` or ``--tau``)
must be provided:

.. code-block:: console

  $ w_timings -W west.h5 -t 100_ps
  Iterations: 50
  Total segments: 9985
  Wall-clock time: 0:00:01.452625
  Simulated physical time ("molecular time"): 5.0 ns
  Aggregate simulation time: 998.5 ns
