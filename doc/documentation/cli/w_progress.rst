w_progress
==========

``w_progress`` shows a live progress dashboard for a WESTPA simulation.

Overview
--------

Usage:

.. code-block:: text

  w_progress [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
             [-W WEST_H5FILE] [--refresh SECONDS]

Tool-specific arguments:

.. code-block:: text

  --refresh SECONDS     Redraw the dashboard every SECONDS seconds (default:
                        1.0).

Run ``w_progress`` in a separate terminal while ``w_run`` runs. It redraws the
dashboard until Ctrl-C.

Examples
--------

.. code-block:: console

  $ w_progress -W west.h5
  WESTPA progress for west.h5 (updated 21:09:25)

  Last written:              2026-06-06 21:09:10
  Current iteration:         23
  Progress:                  22 / 100 iterations (22.0%)
  Segments complete:         41 / 80
  Segments failed:           0
  Recent walltimes:          0:05:22, 0:05:19, 0:05:17, 0:05:21, 0:05:23
  Avg iteration time:        0:05:20
  Completed walltime:        1:57:32
  Completed segments:        1760
  ETA:                       6:56:31

Notes
-----

``w_run`` keeps ``west.h5`` open for writing. HDF5 file locking blocks other
programs from reading it. ``w_progress`` opens the file with locking off, so it
only sees data ``w_run`` has flushed to disk. ``w_run`` flushes at the end of
each iteration and every ``flush_period`` seconds (60 by default) as segments
finish. ``Last written`` shows the time of the last flush.

ETA is iterations left times the average walltime of the last five completed
iterations. ETA shows ``unknown`` until one iteration finishes, or if
``max_total_iterations`` is not set in ``west.cfg``.

If a read fails, for example mid flush, the last good values stay on screen
with the error below. The next refresh tries again.
