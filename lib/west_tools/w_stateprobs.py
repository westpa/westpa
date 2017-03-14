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

from w_direct import DStateProbs
import sys, argparse, os
import work_managers

# Just a shim to make sure everything works and is backwards compatible.
# We're making sure it has the appropriate functions so that it can be called
# as a regular tool, and not a subcommand.

class WStateProbs(DStateProbs):
    subcommand = 'trace'
    help_text = 'averages and CIs for path-tracing kinetics analysis'
    default_output_file = 'stateprobs.h5'
    # This isn't strictly necessary, but for the moment, here it is.
    # We really need to modify the underlying class so that we don't pull this sort of stuff if it isn't necessary.
    # That'll take some case handling, which is fine.
    default_kinetics_file = 'assign.h5'

class WDirect(WESTMasterCommand, WESTParallelTool):
    prog='w_stateprobs'
    subcommands = [WStateProbs]
    subparsers_title = 'calculate state-to-state kinetics by tracing trajectories'
    description = '''\
Calculate average populations and associated errors in state populations from
weighted ensemble data. Bin assignments, including macrostate definitions,
are required. (See "w_assign --help" for more information).

-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, usually "stateprobs.h5") contains the following
dataset:

  /avg_state_pops [state]
    (Structured -- see below) Population of each state across entire
    range specified.

If --evolution-mode is specified, then the following additional dataset is
available:

  /state_pop_evolution [window][state]
    (Structured -- see below). State populations based on windows of
    iterations of varying width.  If --evolution-mode=cumulative, then
    these windows all begin at the iteration specified with
    --start-iter and grow in length by --step-iter for each successive 
    element. If --evolution-mode=blocked, then these windows are all of
    width --step-iter (excluding the last, which may be shorter), the first
    of which begins at iteration --start-iter.
    
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

if __name__ == '__main__':
    print('WARNING: {} is being deprecated.  Please use w_direct instead.'.format(WDirect.prog))
    # If we're not really supporting subcommands...
    import sys
    try:
        if sys.argv[1] != 'trace':
            sys.argv.insert(1, 'trace')
    except:
        sys.argv.insert(1, 'trace')
    WDirect().main()
