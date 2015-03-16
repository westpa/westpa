Analysis
========

Gauging simulation progress and convergence
-------------------------------------------

Progress coordinate distribution (w_pcpdist)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

w_pcpdist and plothist

Kinetics for source/sink simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

w_fluxanl

Kinetics for arbitrary state definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to calculate rate constants, it is necessary to run three different
tools::

- :ref:`w_assign`
- :ref:`w_kinetics`
- :ref:`w_kinavg`

The w_assign tool assigns trajectories to states (states which correspond to a
target bin) at a sub-tau resolution. This allows w_kinetics to properly trace
the trajectories and prepare the data for further analysis.

Although the bin and state definitions can be pulled from the system, it is
frequently more convenient to specify custom bin boundaries and states; this
eliminates the need to know what constitutes a state prior to starting the
simulation. Both files must be in the YAML format, of which there are numerous
examples of online. A quick example for each file follows::

  States:
  ---
  states:
    - label: unbound
      coords:
        - [25,0]
    - label: boun
      coords:
        - [1.5,33.0]

  Bins:
  ---
  bins:
    type: RectilinearBinMapper
    boundaries: [[0.0,1.57,25.0,10000],[0.0,33.0,10000]]

This system has a two dimensional progress coordinate, and two definite states,
as defined by the PMF. The binning used during the simulation was significantly
more complex; defining a smaller progress coordinate (in which we have three
regions: bound, unbound, and in between) is simply a matter of convenience.
Note that these custom bins do not change the simulation in any fashion; you
can adjust state definitions and bin boundaries at will without altering the
way the simulation runs.

The help definition, included by running::
  
  w_assign --help

usually contains the most up-to-date help information, and so more
information about command line options can be obtained from there. To
run with the above YAML files, assuming they are named STATES and BINS,
you would run the following command::

  w_assign --states-from-file STATES --bins-from-file BINS

By default, this produces a .h5 file (named assign.h5); this can be changed via
the command line.

The w_kinetics tool uses the information generated from w_assign to trace
through trajectories and calculate flux with included color information. There
are two main methods to run w_kinetics::

  w_kinetics trace
  w_kinetics matrix

The matrix method is still in development; at this time, trace is the
recommended method.

Once the w_kinetics analysis is complete, you can check for convergence of the
rate constants. WESTPA includes two tools to help you do this: w_kinavg and
ploterr. First, begin by running the following command (keep in mind that
w_kinavg has the same type of analysis as w_kinetics does; whatever method you
chose (trace or matrix) in the w_kinetics step should be used here, as well)::

  w_kinavg trace -e cumulative

This instructs w_kinavg to produce a .h5 file with the cumulative rate
information; by then using ploterr, you can determine whether the rates
have stopped changing::

  ploterr kinavg

By default, this produces a set of .pdf files, containing cumulative rate and
flux information for each state-to-state transition as a function of the WESTPA
iteration. Determine at which iteration the rate stops changing; then, rerun
w_kinavg with the following systems::

  w_kinavg trace --first-iter ITER

where ITER is the beginning of the unchanging region. This will then
output information much like the following::

  fluxes into macrostates:
  unbound: mean=1.712580005863456e-02 CI=(1.596595628304422e-02, 1.808249529394858e-02) * tau^-1
  bound  : mean=5.944989301935855e-04 CI=(4.153556214886056e-04, 7.789568983584020e-04) * tau^-1
  
  fluxes from state to state:
  unbound -> bound  : mean=5.944989301935855e-04 CI=(4.253003401668849e-04, 7.720997503648696e-04) * tau^-1
  bound   -> unbound: mean=1.712580005863456e-02 CI=(1.590547796439216e-02, 1.808154616175579e-02) * tau^-1
  
  rates from state to state:
  unbound -> bound  : mean=9.972502012305491e-03 CI=(7.165030136921814e-03, 1.313767180582492e-02) * tau^-1
  bound   -> unbound: mean=1.819520888349874e-02 CI=(1.704608273094848e-02, 1.926165865735958e-02) * tau^-1

Divide by tau to calculate your rate constant.
