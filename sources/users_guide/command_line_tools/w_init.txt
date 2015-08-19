.. _w_init:

w_init
======

``w_init`` initializes the weighted ensemble simulation, creates the
main HDF5 file and prepares the first iteration.

Overview
--------

Usage::

  $WEST_ROOT/bin/w_init [-h] [-r RCFILE] [--quiet | --verbose | --debug] [--version]
               [--force] [--bstate-file BSTATE_FILE] [--bstate BSTATES]
               [--tstate-file TSTATE_FILE] [--tstate TSTATES]
               [--segs-per-state N] [--no-we] [--wm-work-manager WORK_MANAGER]
               [--wm-n-workers N_WORKERS] [--wm-zmq-mode MODE]
               [--wm-zmq-info INFO_FILE] [--wm-zmq-task-endpoint TASK_ENDPOINT]
               [--wm-zmq-result-endpoint RESULT_ENDPOINT]
               [--wm-zmq-announce-endpoint ANNOUNCE_ENDPOINT]
               [--wm-zmq-heartbeat-interval INTERVAL]
               [--wm-zmq-task-timeout TIMEOUT] [--wm-zmq-client-comm-mode MODE]

Initialize a new WEST simulation, creating the WEST HDF5 file and preparing the
first iteration's segments. Initial states are generated from one or more
"basis states" which are specified either in a file specified with
``--bstates-from``, or by one or more ``--bstate`` arguments. If neither
``--bstates-from`` nor at least one ``--bstate`` argument is provided, then a
default basis state of probability one identified by the state ID zero and
label "basis" will be created (a warning will be printed in this case, to
remind you of this behavior, in case it is not what you wanted). Target states
for (non- equilibrium) steady-state simulations are specified either in a file
specified with ``--tstates-from``, or by one or more ``--tstate`` arguments. If
neither ``--tstates-from`` nor at least one ``--tstate`` argument is provided,
then an equilibrium simulation (without any sinks) will be performed.

Command-Line Options
--------------------

See the `general command-line tool reference <UserGuide:ToolRefs>`_ for more
information on the general options.

State Options
~~~~~~~~~~~~~

::

  --force
    Overwrites any existing simulation data

  --bstate BSTATES
    Add the given basis state (specified as a string
    'label,probability[,auxref]') to the list of basis states (after
    those specified in --bstates-from, if any). This argument may be
    specified more than once, in which case the given states are
    appended in the order they are given on the command line.

  --bstate-file BSTATE_FILE, --bstates-from BSTATE_FILE
    Read basis state names, probabilities, and (optionally) data
    references from BSTATE_FILE.

  --tstate TSTATES
    Add the given target state (specified as a string
    'label,pcoord0[,pcoord1[,...]]') to the list of target states (after
    those specified in the file given by --tstates-from, if any). This
    argument may be specified more than once, in which case the given
    states are appended in the order they appear on the command line.

  --tstate-file TSTATE_FILE, --tstates-from TSTATE_FILE
    Read target state names and representative progress coordinates from
    TSTATE_FILE.

  --segs-per-state N
    Initialize N segments from each basis state (default: 1).

  --no-we, --shotgun
    Do not run the weighted ensemble bin/split/merge algorithm on
    newly-created segments.

Examples
--------

(TODO: write 3 examples; Setting up the basis states, explanation of
bstates and istates. Setting up an equilibrium simulation, w/o target(s)
for recycling. Setting up a simulation with one/multiple target states.)
