.. _running:

Running
=======

Overview
--------

The **w_run** command is used to run weighted ensemble simulations
`configured <setup>` with **w_init**.

Setting simulation limits
-------------------------

Running a simulation
--------------------

Running on a single node
~~~~~~~~~~~~~~~~~~~~~~~~

Running on multiple nodes with MPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running on multiple nodes with ZeroMQ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Managing data
-------------

Recovering from errors
----------------------

By default, information about simulation progress is stored in
**west-JOBID.log** (where JOBID refers to the job ID given by the submission
engine); any errors will be logged here.

The :ref:`w_progress` command can also be used to monitor a running or recently
updated simulation. It reports the current iteration, segment
prepared and failed counts, recent iteration timings, and an estimated time to
completion when available. During active runs, ``w_progress`` reads the
``west.h5.progress.json`` status file written by ``w_run``

- The error "could not read pcoord from 'tempfile': progress coordinate has
  incorrect shape" may come about from multiple causes; it is possible that the
  progress coordinate length is incorrectly specified in system.py
  (**self.pcoord_len**), or that GROMACS (or whatever simulation package you
  are using) had an error during the simulation.

- The first case will be obvious by what comes after the message: (XX, YY)
  (where XX is non-zero), expected (ZZ, GG) (whatever is in system.py). This
  can be corrected by adjusting system.py.
- In the second case, the progress coordinate length is 0; this
  indicates that no progress coordinate data exists (null string), which
  implies that the simulation software did not complete successfully. By
  default, the simulation package (GROMACS or otherwise) terminal output is
  stored in a log file inside of seg_logs. Any error that occurred during the
  actual simulation will be logged here, and can be corrected as needed.
