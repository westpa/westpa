.. _w_stateprobs:

w_stateprobs
============

WARNING: w_stateprobs is being deprecated.  Please use w_direct instead.

usage::

 w_stateprobs trace [-h] [-W WEST_H5FILE] [--first-iter N_ITER] [--last-iter N_ITER]
                          [--step-iter STEP] [-a ASSIGNMENTS] [-o OUTPUT] [-k KINETICS]
                          [--disable-bootstrap] [--disable-correl] [--alpha ALPHA]
                          [--autocorrel-alpha ACALPHA] [--nsets NSETS] [-e {cumulative,blocked,none}]
                          [--window-frac WINDOW_FRAC] [--disable-averages]

Calculate average populations and associated errors in state populations from
weighted ensemble data. Bin assignments, including macrostate definitions,
are required. (See "w_assign --help" for more information).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "direct.h5") contains the following
dataset::

  /avg_state_probs [state]
    (Structured -- see below) Population of each state across entire
    range specified.

  /avg_color_probs [state]
    (Structured -- see below) Population of each ensemble across entire
    range specified.

If --evolution-mode is specified, then the following additional datasets are
available::

  /state_pop_evolution [window][state]
    (Structured -- see below). State populations based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.

  /color_prob_evolution [window][state]
    (Structured -- see below). Ensemble populations based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.

The structure of these datasets is as follows::

  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).

  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).

  expected
    (Floating-point) Expected (mean) value of the observable as evaluated within
    this window, in units of inverse tau.

  ci_lbound
    (Floating-point) Lower bound of the confidence interval of the observable
    within this window, in units of inverse tau.

  ci_ubound
    (Floating-point) Upper bound of the confidence interval of the observable
    within this window, in units of inverse tau.

  stderr
    (Floating-point) The standard error of the mean of the observable
    within this window, in units of inverse tau.

  corr_len
    (Integer) Correlation length of the observable within this window, in units
    of tau.

Each of these datasets is also stamped with a number of attributes::

  mcbs_alpha
    (Floating-point) Alpha value of confidence intervals. (For example,
    *alpha=0.05* corresponds to a 95% confidence interval.)

  mcbs_nsets
    (Integer) Number of bootstrap data sets used in generating confidence
    intervals.

  mcbs_acalpha
    (Floating-point) Alpha value for determining correlation lengths.

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
                        Store results in OUTPUT (default: stateprobs.h5).

input/output options::

  -k KINETICS, --kinetics KINETICS
                        Populations and transition rates are stored in KINETICS (default: assign.h5).

confidence interval calculation options::

  --disable-bootstrap, -db
                        Enable the use of Monte Carlo Block Bootstrapping.
  --disable-correl, -dc
                        Disable the correlation analysis.
  --alpha ALPHA         Calculate a (1-ALPHA) confidence interval' (default: 0.05)
  --autocorrel-alpha ACALPHA
                        Evaluate autocorrelation to (1-ACALPHA) significance. Note that too small an
                        ACALPHA will result in failure to detect autocorrelation in a noisy flux signal.
                        (Default: same as ALPHA.)
  --nsets NSETS         Use NSETS samples for bootstrapping (default: chosen based on ALPHA)

calculation options::

  -e {cumulative,blocked,none}, --evolution-mode {cumulative,blocked,none}
                        How to calculate time evolution of rate estimates. ``cumulative`` evaluates rates
                        over windows starting with --start-iter and getting progressively wider to --stop-
                        iter by steps of --step-iter. ``blocked`` evaluates rates over windows of width
                        --step-iter, the first of which begins at --start-iter. ``none`` (the default)
                        disables calculation of the time evolution of rate estimates.
  --window-frac WINDOW_FRAC
                        Fraction of iterations to use in each window when running in ``cumulative`` mode.
                        The (1 - frac) fraction of iterations will be discarded from the start of each
                        window.

misc options::

  --disable-averages, -da
                        Whether or not the averages should be printed to the console (set to FALSE if flag
                        is used).