from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math
import numpy, h5py
import wemd, wemdtools
from wemdtools.transitions.transacc import TransitionEventAccumulator
from wemdtools.files import load_npy_or_text

import logging
log = logging.getLogger('w_ttimes')


def pstatus(*pargs, **kwargs):
    global args
    if not args.quiet_mode:
        print(*pargs, file=sys.stderr, **kwargs)

#parser = wemd.rc.common_arg_parser('w_ttimes', description='''
#Perform lifetime, transition, and kinetic analysis on brute force data.''')

parser = argparse.ArgumentParser()
parser.add_argument('input_files', nargs='+', metavar='INPUT_FILE',
                    help='Input data. Specify one or more text files or numpy (.npy) files.')
parser.add_argument('--usecols', metavar='COLS',
                    help='Use only COLS (a comma-separated list of zero-based column indices) from the given input.')
#parser.add_argument('--cat', action='store_true', help='Multiple input files are parts of a single trajectory')
parser.add_argument('-H', '--h5', metavar='H5FILE:H5GROUP', default='analysis.h5:/w_ttimes',
                    help='Store intermediate values and results in H5FILE, in group H5GROUP. '
                        +'(Default: analysis.h5:/w_ttimes)')
parser.add_argument('--record-all-crossings', action='store_true',
                    help='Record all boundary crossing events in the HDF5 file. This dramatically increases processing '
                        +'time and disk space. Boundary crossing counts are calculated regardless of this option.')
parser.add_argument('--record-self-transitions', action='store_true',
                    help='Record all self-transition events in the HDF5 file. This dramatically increases processing '
                        +'time and disk space. Self-transition event counts are calculated regardless of this option.')
parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit status messages')
wemdtools.bins.add_region_set_options(parser)
args = parser.parse_args()

region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stderr)
n_bins = len(region_set.get_all_bins())

try:
    h5fn, h5gn = args.h5.rsplit(':',1)
except ValueError:
    # No colon
    h5fn = args.h5
    h5gn = 'w_ttimes'
    
if args.usecols:
    args.usecols = list(map(int,args.usecols.split(',')))
    pstatus('Using columns', args.usecols)

h5file = h5py.File(h5fn)
try:
    del h5file[h5gn]
except KeyError:
    pass
h5group = h5file.create_group(h5gn)

transacc = TransitionEventAccumulator(n_bins, h5group)
transacc.record_all_crossings = bool(args.record_all_crossings)
transacc.record_self_transitions = bool(args.record_self_transitions)

for input_file in args.input_files:
    pstatus('Processing', input_file)
    
    pstatus('  Loading file...')
    pcoords = load_npy_or_text(input_file)
    if args.usecols:
        pcoords = pcoords[:, args.usecols]
    #h5group['pcoords'] = pcoords
    
    pstatus('  Assigning to bins...')
    assignments = numpy.empty((len(pcoords),), dtype=numpy.min_scalar_type(n_bins))
    assignments[:] = region_set.map_to_all_indices(pcoords)
    #h5group['assignments'] = assignments
    del pcoords
    
    pstatus('  Finding transitions...')
    transacc.start_accumulation(assignments, numpy.ones((len(assignments),)), numpy.ones((len(assignments),n_bins)))
    transacc.flush_transition_data()
    del assignments
    
h5group['n_trans'] = transacc.n_trans
