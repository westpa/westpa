``w_timings``
=============

``w_timings`` prints timing information for a WESTPA simulation.

Overview
--------

Usage:

.. code-block:: text

  w_timings [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
            [-W WEST_H5FILE] [--first-iter N_ITER] [--last-iter N_ITER] [-t TAU]

Tool-specific arguments:

.. code-block:: text

  -t TAU, --tau TAU     WE resampling interval (format: <value>_<unit>, where
                        <value> is a positive integer and <unit> is 'as', 'fs',
                        'ps', 'ns', 'us', 'ms', 's', 'm', 'h', 'D', or 'W').

Examples
--------

The default output includes the wall-clock time and, if per-segment CPU times
were recorded by the propagator, the total CPU time:

.. code-block:: console

  $ w_timings -W west.h5
  Iterations: 100
  Total segments: 7775
  Wall-clock time: 2 days, 4:35:11.749094
  Total CPU time: 17 days, 9:09:00.556985

To determine the maximum trajectory length and aggregate simulation time,
the resampling interval (``TAU``) must be specified:

.. code-block:: console

  $ w_timings -W west.h5 -t 50_ps
  Iterations: 100
  Total segments: 7775
  Wall-clock time: 2 days, 4:35:11.749094
  Total CPU time: 17 days, 9:09:00.556985
  Maximum trajectory length: 5.0 ns
  Aggregate simulation time: 388.75 ns
