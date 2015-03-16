.. _w_run:

w_run
=====

``w_run`` starts or continues a weighted ensemble simualtion.

Overview
--------

Usage::

  $WEST_ROOT/bin/w_run [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
               [--oneseg ] [--wm-work-manager WORK_MANAGER]
               [--wm-n-workers N_WORKERS] [--wm-zmq-mode MODE]
               [--wm-zmq-info INFO_FILE] [--wm-zmq-task-endpoint TASK_ENDPOINT]
               [--wm-zmq-result-endpoint RESULT_ENDPOINT]
               [--wm-zmq-announce-endpoint ANNOUNCE_ENDPOINT]
               [--wm-zmq-heartbeat-interval INTERVAL]
               [--wm-zmq-task-timeout TIMEOUT] [--wm-zmq-client-comm-mode MODE]

Command-Line Options
--------------------

See the :ref:`command-line tool index <command_line_tool_index>` for
more information on the general options.

Segment Options
~~~~~~~~~~~~~~~

::
  --oneseg
    Only propagate one segment (useful for debugging propagators)

Example
-------

A simple example for using w_run (mostly taken from odld example that
is available in the main WESTPA distribution)::

  $WEST_ROOT/bin/w_run &> west.log

This commands starts up a serial weighted ensemble run and pipes the results
into the west.log file. As a side note ``--debug`` option is very useful for
debugging the code if something goes wrong.
