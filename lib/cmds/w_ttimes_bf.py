from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, warnings
import numpy, h5py
import wemd, wemdtools
from wemdtools.transitions.transacc import TransitionEventAccumulator
from wemdtools.files import load_npy_or_text

import logging
log = logging.getLogger('w_ttimes')

ciinfo_dtype = numpy.dtype([('expectation', numpy.float64),
                            ('ci_lower', numpy.float64),
                            ('ci_upper', numpy.float64),
                            ])

def pstatus(*pargs, **kwargs):
    global args
    if not args.quiet_mode:
        print(*pargs, file=sys.stderr, **kwargs)
        
#parser = wemd.rc.common_arg_parser('w_ttimes', description='''
#Perform lifetime, transition, and kinetic analysis on brute force data.''')

parser = argparse.ArgumentParser()
parser.add_argument('input_files', nargs='*', metavar='INPUT_FILE',
                    help='Input data. Specify one or more text files or numpy (.npy) files.')
parser.add_argument('--usecols', metavar='COLS',
                    help='Use only COLS (a comma-separated list of zero-based column indices) from the given input.')
parser.add_argument('--chunksize', type=long, default=100000,
                    help='Process each input in chunks of not more than CHUNKSIZE data points (default: %(default)s). '
                        +'This will only decrease memory use when using numpy (.npy) input files.')
#parser.add_argument('--cat', action='store_true', help='Multiple input files are parts of a single trajectory')
parser.add_argument('-H', '--h5', metavar='H5FILE:H5GROUP', default='analysis.h5:/w_ttimes',
                    help='Store intermediate values and results in H5FILE, in group H5GROUP. '
                        +'(Default: analysis.h5:/w_ttimes)')
parser.add_argument('--statistics-only', action='store_true',
                    help='Do not determine transitions; only process existing data in the HDF5 file.')
parser.add_argument('--record-all-crossings', action='store_true',
                    help='Record all boundary crossing events in the HDF5 file. This dramatically increases processing '
                        +'time and disk space. Boundary crossing counts are calculated regardless of this option.')
parser.add_argument('--record-self-transitions', action='store_true',
                    help='Record all self-transition events (i->j->...->i) in the HDF5 file. This dramatically increases processing '
                        +'time and disk space. Self-transition event counts are calculated regardless of this option.')
parser.add_argument('--dt', type=float, default=1.0,
                    help='Assume input data has a time spacing of DT (default: 1.0)')

parser.add_argument('-E', '--edstats', dest='ed_stats', default='edstats.txt',
                    help='Store event duration statistics in ED_STATS (default: edstats.txt)')
parser.add_argument('-F', '--fptstats', dest='fpt_stats', default='fptstats.txt',
                    help='Store first passage time statistics in FPT_STATS (default: fptstats.txt)')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not include headers in text output files (default: include headers)')
parser.add_argument('-l', '--labels', dest='print_labels', action='store_true',
                    help='Print bin labels in output files, if available (default: do not print bin labels)')

wemdtools.stats.mcbs.add_mcbs_options(parser)

parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='Do not emit status messages')
wemdtools.bins.add_region_set_options(parser)
args = parser.parse_args()


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

# Trace trajectories and record transitions
if not args.statistics_only:
    try:
        del h5file[h5gn]
    except KeyError:
        pass
    h5group = h5file.create_group(h5gn)
        
    region_set = wemdtools.bins.get_region_set_from_args(args, status_stream=sys.stderr)
    n_bins = len(region_set.get_all_bins())
    
    transacc = TransitionEventAccumulator(n_bins, h5group)
    transacc.record_all_crossings = bool(args.record_all_crossings)
    transacc.record_self_transitions = bool(args.record_self_transitions)
        
    for input_file in args.input_files:
        pstatus('Processing', input_file)
        
        pstatus('  Loading file...')
        pcoords = load_npy_or_text(input_file)
        
        for ii in xrange(0, len(pcoords), args.chunksize):
            pstatus('\r  Processing {:d}/{:d}'.format(ii,len(pcoords)), end='')
            
            if args.usecols:
                pcoord_chunk = pcoords[ii:ii+args.chunksize, args.usecols]
            else:
                pcoord_chunk = pcoords[ii:ii+args.chunksize]
            
            assignments = numpy.empty((len(pcoord_chunk),), dtype=numpy.min_scalar_type(n_bins))
            assignments[:] = region_set.map_to_all_indices(pcoord_chunk)
            weights = numpy.ones((len(assignments),))
            binpops = numpy.ones((len(assignments),n_bins))
            
            if ii == 0:
                transacc.start_accumulation(assignments, weights, binpops)
            else:
                transacc.continue_accumulation(assignments, weights, binpops)
            
            transacc.flush_transition_data()
            del pcoord_chunk, assignments
        pstatus('')
        
    h5group['n_trans'] = transacc.n_trans
else:
    h5group = h5file[h5gn]
    region_set = None
    n_bins = h5group['n_trans'].shape[0]

print('Number of transitions:')
print(h5group['n_trans'][...])
    
# Obtain statistics
pstatus('Calculating statistics...')

alpha = 1.0 - args.confidence
n_sets = args.bssize or wemdtools.stats.mcbs.get_bssize(alpha)
syn_avg_durations = numpy.empty((n_sets,), numpy.float64)
syn_avg_fpts = numpy.empty((n_sets,), numpy.float64)

lbi = int(math.floor(n_sets*alpha/2))
ubi = int(math.ceil(n_sets*(1-alpha/2))) 

pstatus('  Using bootstrap sampling of {:d} sets'.format(n_sets))
pstatus('  {:g}% confidence interval (alpha={:g}) requested [bounding indices ({:d},{:d})]'.format(args.confidence*100,
                                                                                                  alpha,
                                                                                                  lbi,
                                                                                                  ubi))

transdat_ds = h5group['transitions']
transdat_ibin = transdat_ds['initial_bin'][:]
transdat_fbin = transdat_ds['final_bin'][:]

durations = numpy.zeros((n_bins,n_bins), ciinfo_dtype)
fpts = numpy.zeros((n_bins,n_bins), ciinfo_dtype)

for ibin in xrange(n_bins):
    trans_ibin = transdat_ds[transdat_ibin == ibin]
    for fbin in xrange(n_bins):
        pstatus('\r  Bin {:d}->{:d}'.format(ibin,fbin), end='')
        trans_ifbins = trans_ibin[trans_ibin['final_bin'] == fbin]

        # Suppress divide-by-zero errors; we actually want the NaNs        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            
            durations[ibin,fbin]['expectation'] = trans_ifbins['duration'].mean() * args.dt
            fpts[ibin,fbin]['expectation'] =      trans_ifbins['fpt'].mean() * args.dt
        
        dlen = len(trans_ifbins)
        if not dlen:
            continue # with next fbin
        
        for iset in xrange(n_sets):
            indices = numpy.random.randint(dlen, size=(dlen,))
            syn_trans = trans_ifbins[indices]
            with_fpts = syn_trans[syn_trans['fpt'] > 0]
            
            syn_avg_durations[iset] = syn_trans['duration'].mean() * args.dt
            syn_avg_fpts[iset] = with_fpts['fpt'].mean() * args.dt
            
        syn_avg_durations.sort()
        syn_avg_fpts.sort()
        
        durations[ibin,fbin]['ci_lower'] = syn_avg_durations[lbi]
        durations[ibin,fbin]['ci_upper'] = syn_avg_durations[ubi]
        
        fpts[ibin,fbin]['ci_lower'] = syn_avg_fpts[lbi]
        fpts[ibin,fbin]['ci_upper'] = syn_avg_fpts[ubi]
pstatus('')
for dsname in ('duration', 'fpt'):
    try:
        del h5group[dsname]
    except KeyError:
        pass

h5group['duration'] = durations
h5group['fpt'] = fpts

for dataset in (h5group['duration'], h5group['fpt']):
    dataset.attrs['dt'] = args.dt
    dataset.attrs['alpha'] = alpha

max_ibin_width = len(str(n_bins))-1
format_2d = '{ibin:{mw}d}    {fbin:{mw}d}    {0:20.15g}    {1:20.15g}    {2:20.15g}    {3:20.15g}    {4:20.15g}    {5:20.15g}\n'
for (filename, data, title) in ((args.ed_stats, durations, 'event duration'),
                                (args.fpt_stats, fpts, 'first passage time')):
    if not filename: continue
    outfile = open(filename, 'wt')
    if not args.suppress_headers:
        outfile.write('''\
# {title:} statistics
# confidence interval = {confidence:g}%
# ----
# column 0: initial bin index
# column 1: final bin index
# column 2: lower bound of confidence interval
# column 3: upper bound of confidence interval
# column 4: width of confidence interval
# column 5: relative width of confidence interval [abs(width/average)]
# column 6: symmetrized error [max(upper-average, average-lower)]
# ----
'''.format(title=title, confidence=args.confidence*100))
        if region_set is not None and args.print_labels:
            wemdtools.bins.print_labels(region_set, outfile)
            outfile.write('----\n')
    
    for ibin in xrange(n_bins):
        for fbin in xrange(n_bins):
            mean = data[ibin,fbin]['expectation']
            lb = data[ibin,fbin]['ci_lower']
            ub = data[ibin,fbin]['ci_upper']
            ciwidth = ub - lb
            relciwidth = abs(ciwidth/mean)
            symmerr = max(mean-lb,ub-mean)
            
            outfile.write(format_2d.format(*map(float,(mean,lb,ub,ciwidth,relciwidth,symmerr)),
                                           ibin=ibin,fbin=fbin,mw=max_ibin_width))

