from westpa.tools import WESTMasterCommand, WESTParallelTool
from warnings import warn

from westpa.cli.tools.w_reweight import RWAverage

# Just a shim to make sure everything works and is backwards compatible.
# We're making sure it has the appropriate functions so that it can be called
# as a regular tool, and not a subcommand.


class PAAverage(RWAverage):
    subcommand = 'average'
    help_text = ''
    default_output_file = 'kinrw.h5'
    # This isn't strictly necessary, but for the moment, here it is.
    # We really need to modify the underlying class so that we don't pull this sort of stuff if it isn't necessary.
    # That'll take some case handling, which is fine.
    default_kinetics_file = 'flux_matrices.h5'


class WReweight(WESTMasterCommand, WESTParallelTool):
    prog = 'w_postanalysis_reweight'
    subcommands = [PAAverage]
    subparsers_title = 'calculate state-to-state kinetics by tracing trajectories'
    description = '''\
A convenience function to run kinetics/probs. Bin assignments,
including macrostate definitions, are required. (See
"w_assign --help" for more information).

For more information on the individual subcommands this subs in for, run
w_reweight {kinetics/probs} --help.

-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''


def entry_point():
    warn('{} is being deprecated.  Please use w_reweight instead.'.format(WReweight.prog))
    # If we're not really supporting subcommands...
    import sys

    try:
        if sys.argv[1] != 'average':
            sys.argv.insert(1, 'average')
    except Exception:
        sys.argv.insert(1, 'average')
    WReweight().main()


if __name__ == '__main__':
    entry_point()
