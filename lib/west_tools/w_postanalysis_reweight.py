# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from westtools import (WESTMasterCommand, WESTParallelTool, WESTDataReader, IterRangeSelection, WESTSubcommand, WESTToolComponent, WESTTool,
                       ProgressIndicatorComponent)

from w_reweight import RWAverage
import sys, argparse, os
import work_managers

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
    prog='w_postanalysis_reweight'
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

if __name__ == '__main__':
    print('WARNING: {} is being deprecated.  Please use w_reweight instead.'.format(WReweight.prog))
    # If we're not really supporting subcommands...
    import sys
    try:
        if sys.argv[1] != 'average':
            sys.argv.insert(1, 'average')
    except:
        sys.argv.insert(1, 'average')
    WReweight().main()
