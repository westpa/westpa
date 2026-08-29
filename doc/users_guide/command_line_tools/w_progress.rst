.. _w_progress:

w_progress
==========

``w_progress`` shows a live terminal dashboard for a running or recently updated
WESTPA simulation.

Usage:

.. code-block:: text

  w_progress [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
             [-W WEST_H5FILE] [--refresh SECONDS]

Examples
--------

Refresh the dashboard once per second using the default WEST data file:

.. code-block:: console

  $ w_progress

Read a specific WEST HDF5 file:

.. code-block:: console

  $ w_progress -W west.h5

Refresh every 10 seconds:

.. code-block:: console

  $ w_progress --refresh 10

Output
------

The dashboard reports the current iteration, latest completed iteration,
requested total iterations when available, prepared and failed segment counts
for the current iteration with the current segment total, recent iteration
wall-clock times, completed wall-clock time when available, completed segment
count when available, and a simple ETA based on the last five completed
iterations.

``w_progress`` is read-only. It does not control or modify a WESTPA simulation.
The command is intended to be run in a separate terminal while a simulation is
running. By default, it clears and redraws the terminal on each refresh.

During active runs started with the current ``w_run``, live values come from a
small sidecar status file named ``west.h5.progress.json``. When that file is not
available, or when the run has completed, ``w_progress`` falls back to reading
``west.h5`` directly.

The top ``Updated`` timestamp shows when ``w_progress`` refreshed the dashboard.
During active runs, ``Live status updated`` shows when ``w_run`` last wrote new
live status data.
