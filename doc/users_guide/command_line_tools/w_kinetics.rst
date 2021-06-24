.. _w_kinetics:

w_kinetics
==========

WARNING: w_kinetics is being deprecated.  Please use w_direct instead.

usage::

 w_kinetics trace [-h] [-W WEST_H5FILE] [--first-iter N_ITER] [--last-iter N_ITER]
                        [--step-iter STEP] [-a ASSIGNMENTS] [-o OUTPUT]

Calculate state-to-state rates and transition event durations by tracing
trajectories.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).

This subcommand for w_direct is used as input for all other w_direct
subcommands, which will convert the flux data in the output file into
average rates/fluxes/populations with confidence intervals.

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "direct.h5") contains the
following datasets::

  ``/conditional_fluxes`` [iteration][state][state]
    *(Floating-point)* Macrostate-to-macrostate fluxes. These are **not**
    normalized by the population of the initial macrostate.

  ``/conditional_arrivals`` [iteration][stateA][stateB]
    *(Integer)* Number of trajectories arriving at state *stateB* in a given
    iteration, given that they departed from *stateA*.

  ``/total_fluxes`` [iteration][state]
    *(Floating-point)* Total flux into a given macrostate.

  ``/arrivals`` [iteration][state]
    *(Integer)* Number of trajectories arriving at a given state in a given
    iteration, regardless of where they originated.

  ``/duration_count`` [iteration]
    *(Integer)* The number of event durations recorded in each iteration.

  ``/durations`` [iteration][event duration]
    *(Structured -- see below)*  Event durations for transition events ending
    during a given iteration. These are stored as follows:

      istate
        *(Integer)* Initial state of transition event.
      fstate
        *(Integer)* Final state of transition event.
      duration
        *(Floating-point)* Duration of transition, in units of tau.
      weight
        *(Floating-point)* Weight of trajectory at end of transition, **not**
        normalized by initial state population.

Because state-to-state fluxes stored in this file are not normalized by
initial macrostate population, they cannot be used as rates without further
processing. The ``w_direct kinetics`` command is used to perform this normalization
while taking statistical fluctuation and correlation into account. See
``w_direct kinetics --help`` for more information.  Target fluxes (total flux
into a given state) require no such normalization.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------

optional arguments::

  -h, --help            show this help message and exit

WEST input data options::

  -W WEST_H5FILE, --west-data WEST_H5FILE
                        Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in
                        west.cfg).

iteration range::

  --first-iter N_ITER   Begin analysis at iteration N_ITER (default: 1).
  --last-iter N_ITER    Conclude analysis with N_ITER, inclusive (default: last completed iteration).
  --step-iter STEP      Analyze/report in blocks of STEP iterations.

input/output options::

  -a ASSIGNMENTS, --assignments ASSIGNMENTS
                        Bin assignments and macrostate definitions are in ASSIGNMENTS (default:
                        assign.h5).
  -o OUTPUT, --output OUTPUT
                        Store results in OUTPUT (default: kintrace.h5).