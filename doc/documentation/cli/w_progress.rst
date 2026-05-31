w_progress
==========

``w_progress`` shows a live terminal dashboard for a WESTPA simulation.

Overview
--------

Usage:

.. code-block:: text

  w_progress [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
             [-W WEST_H5FILE] [--refresh SECONDS]

Tool-specific arguments:

.. code-block:: text

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE.
  --refresh SECONDS     Refresh the dashboard every SECONDS seconds.

Examples
--------

Monitor the default ``west.h5`` file, refreshing once per second:

.. code-block:: console

  $ w_progress

Monitor a specific WEST HDF5 file:

.. code-block:: console

  $ w_progress -W /path/to/west.h5

Refresh every 10 seconds:

.. code-block:: console

  $ w_progress --refresh 10

Notes
-----

``w_progress`` reads data in this order:

#. If the sidecar file written by ``w_run`` exists and reports an active run
   state, render that live status without opening ``west.h5``.
#. Otherwise, read ``west.h5`` directly for the most complete stopped-run or
   completed-run summary.
#. If ``west.h5`` cannot be read and the sidecar reports a completed run, render
   the completed sidecar status as a fallback.

The sidecar is written next to the WEST HDF5 file, for example
``west.h5.progress.json``. This avoids opening ``west.h5`` while the running
simulation owns the HDF5 file lock.

The top ``Updated`` timestamp is the time when ``w_progress`` refreshed the
dashboard. For active runs, the ``Live status updated`` row shows when ``w_run``
last wrote the sidecar data.

``w_progress`` is intended to be run in a separate terminal while a simulation
is running. By default, it clears and redraws the terminal on each refresh.

When the live status file is unavailable, or when the run has completed,
``w_progress`` falls back to reading the latest available state from
``west.h5``.

The ETA is a simple estimate based on the requested total iterations from
``west.cfg`` and the average wall-clock time of recent completed iterations. If
either value is unavailable, the ETA is shown as ``unknown``.
