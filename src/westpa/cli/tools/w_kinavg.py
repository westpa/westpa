from westpa.tools import WESTMasterCommand, WESTParallelTool

from westpa.cli.tools.w_direct import DKinAvg

# Just a shim to make sure everything works and is backwards compatible.


class WKinAvg(DKinAvg):
    subcommand = 'trace'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_kinetics_file = 'kintrace.h5'
    default_output_file = 'kinavg.h5'


class WDirect(WESTMasterCommand, WESTParallelTool):
    prog = 'w_kinavg'
    subcommands = [WKinAvg]
    subparsers_title = 'direct kinetics analysis schemes'
    description = '''\
Calculate average rates and associated errors from weighted ensemble data. Bin
assignments (usually "assignments.h5") and kinetics data (usually
"kintrace.h5" or "kinmat.h5") data files must have been previously generated
(see "w_assign --help" and "w_kinetics --help" for information on generating
these files).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "kinavg.h5") contains the following
dataset:

  /avg_rates [state,state]
    (Structured -- see below) State-to-state rates based on entire window of
    iterations selected.

For trace mode, the following additional datasets are generated:

  /avg_total_fluxes [state]
    (Structured -- see below) Total fluxes into each state based on entire
    window of iterations selected.

  /avg_conditional_fluxes [state,state]
    (Structured -- see below) State-to-state fluxes based on entire window of
    iterations selected.

If --evolution-mode is specified, then the following additional dataset is
available:

  /rate_evolution [window][state][state]
    (Structured -- see below). State-to-state rates based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.

If --evolution-mode is specified in trace mode, the following additional
datasets are available:

  /target_flux_evolution [window,state]
    (Structured -- see below). Total flux into a given macro state based on
    windows of iterations of varying width, as in /rate_evolution.

  /conditional_flux_evolution [window,state,state]
    (Structured -- see below). State-to-state fluxes based on windows of
    varying width, as in /rate_evolution.

The structure of these datasets is as follows:

  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).

  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).

  expected
    (Floating-point) Expected (mean) value of the rate as evaluated within
    this window, in units of inverse tau.

  ci_lbound
    (Floating-point) Lower bound of the confidence interval on the rate
    within this window, in units of inverse tau.

  ci_ubound
    (Floating-point) Upper bound of the confidence interval on the rate
    within this window, in units of inverse tau.

  corr_len
    (Integer) Correlation length of the rate within this window, in units
    of tau.

Each of these datasets is also stamped with a number of attributes:

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
'''


def entry_point():
    print('WARNING: {} is being deprecated.  Please use w_direct instead.'.format(WDirect.prog))
    import sys

    try:
        if sys.argv[1] != 'trace':
            sys.argv.insert(1, 'trace')
    except Exception:
        sys.argv.insert(1, 'trace')
    WDirect().main()


if __name__ == '__main__':
    entry_point()
