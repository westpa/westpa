from __future__ import division, print_function
import os, sys
import numpy

import logging
log = logging.getLogger('w_pcpdist')
import wemd, wemdtools

try:
    import matplotlib
except ImportError:
    pyplot = None
else:
    try:
        matplotlib.use('PDF')
        from matplotlib import pyplot
    except Exception as e:
        log.info('could not select matplotlib PDF backend: {}'.format(e))
        pyplot = None
        
    # Color map based on that in Hovmoeller et al. "Conformations of amino acids in proteins." 
    # Act. Cryst. D., 2002.58. doi:10.1107/S0907444902003359
    cmap_data = numpy.array( [ (124,   0,  24),
                               (211,   0,  32),
                               (244, 184,   0),
                               (245, 235,   0),
                               (129, 183,   2),
                               ( 32, 128,  38),
                               ( 21,  27,  87),
                               ( 36,  62, 137),
                               (178, 220, 245),
                               (255, 255, 255) ], numpy.float32)
    cmap_data /= 255.0
    cm_hovmol   = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol', cmap_data)
    cm_hovmol_r = matplotlib.colors.LinearSegmentedColormap.from_list('hovmol_r', numpy.flipud(cmap_data))

        
parser = wemd.rc.common_arg_parser(description = '''\
Calculate the probability distribution of progress coordinates. Separate output is generated
for each dimension of the progress coordinate.  All data is returned as probability density, 
rather than probability (i.e. each resulting histogram integrates to unity, rather than adds
to unity).
''')

# Subset options
parser.add_argument('-b', '--begin', '--start', dest='start_iter', type=int, default=1,
                    help='Begin accumulation at iteration START_ITER (default: first iteration)')
parser.add_argument('-e', '--end', '--stop', dest='stop_iter', type=int,
                    help='Accumulate through iteration STOP_ITER (default: last completed iteration)')

# Histogram options
parser.add_argument('-B', '--bins', dest='nbins', type=int, default=1000,
                    help='Use NBINS bins for each histogram (default: 1000)')

# Output options
parser.add_argument('-o', '--output', dest='output_pattern', default='pcprob_%d.txt',
                    help='Store average distributions in OUTPUT_PATTERN, which must contain a formatting code which will be '
                        +'replaced by the index of the progress coordinate dimension (default: pcprob_%%d.txt).')
parser.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                    help='Do not write headers to output files (default: write headers).')
parser.add_argument('--quiet', dest='quiet_mode', action='store_true',
                    help='''Do not emit periodic status messages (default: emit status messages if standard output
                    is a terminal)''')
parser.add_argument('datafile', nargs='?',
                    help='Read progress coordinate from DATAFILE(s) (default: load WEMD HDF5 file specified in wemd.cfg).')

# Plotting options
if pyplot:
    parser.add_argument('--plot', dest='plot_pattern', default='pcprob_%d.pdf',
                        help='Store a plot of the probability distribution of each progress coordinate dimension '
                            +'in PLOT_PATTERN, which must contain a formatting code which will be '
                            +'replaced by the progress coordinate dimension (default: pcprob_%%d.pdf)')
    parser.add_argument('--timeplot-lin', dest='timeplot_pattern_lin', default='pcprob_evol_%d_lin.pdf',
                        help='Store a colormap plot of the evolution of the probability distribution of '
                            +'each coordinate dimension in TIMEPLOT_PATTERN_LIN, which must contain a formatting ' 
                            +'code which will be replaced by the progress coordinate dimension  '
                            +'(default: pcprob_evol_%%d_lin.pdf)')
    parser.add_argument('--timeplot-log', dest='timeplot_pattern_log', default='pcprob_evol_%d_log.pdf',
                        help='Store a colormap plot of the evolution of the log of the probability distribution of '
                            +'each coordinate dimension in TIMEPLOT_PATTERN_LOG, which must contain a formatting ' 
                            +'code which will be replaced by the progress coordinate dimension  '
                            +'(default: pcprob_evol_%%d_log.pdf)')    
    parser.add_argument('--vrange', dest='linear_vrange',
                        help='Use LINEAR_VRANGE (a comma-separated pair of values) as the range for the color bar '
                            +'for the linear-scale probability evolution plot. (Default: range of probability density.)')
    parser.add_argument('--lvrange', dest='log_vrange',
                        help='Use LOG_VRANGE (a comma-separated pair of values) as the range for the color bar '
                            +'for the log-scale probability evolution plot. (Default: range of log of nonzero probability density)')    

args = parser.parse_args()

if pyplot:
    if args.linear_vrange:
        vmin,vmax = map(float,args.linear_vrange.split(','))
    if args.log_vrange:
        lvmin,lvmax = map(float,args.log_vrange.split(','))
            
wemd.rc.config_logging(args, 'w_pcpdist')
runtime_config = wemd.rc.read_config(args.run_config_file)
runtime_config.update_from_object(args)
if args.datafile:
    from wemd.util.config_dict import ConfigDict
    if runtime_config is None: runtime_config = ConfigDict()
    runtime_config['data.h5file'] = args.datafile
    print('reading WEMD data from', args.datafile)
else:
    print('reading data from WEMD simulation')
    if runtime_config is None:
        runtime_config = wemd.rc.read_config(args.run_config_file)
sim_manager = wemd.rc.load_sim_manager(runtime_config)
sim_manager.load_data_manager()
data_manager = wemdtools.data_manager.CachingDataReader(sim_manager.data_manager)

start_iter = args.start_iter
stop_iter = data_manager.current_iteration - 1
if args.stop_iter and args.stop_iter <= stop_iter:
    stop_iter = args.stop_iter
n_iters = stop_iter - start_iter + 1

pcoord_dtype = data_manager.get_pcoord_dataset(start_iter).dtype
# 1 dim of seg_id, 1 dim of time, remaining dims are pcoord
pcoord_ndim  = len(data_manager.get_pcoord_dataset(start_iter).shape)-2

# For each dimension...
#    Find min and max pcoord value
#    Choose a resolution
#    Bin data
#    Do weighted histogram

min_pcoord = numpy.empty((pcoord_ndim,), dtype=pcoord_dtype)
max_pcoord = numpy.empty((pcoord_ndim,), dtype=pcoord_dtype)
#var_pcoord = numpy.empty((pcoord_ndim,), dtype=numpy.float64)
n_pcoord_points = numpy.empty((pcoord_ndim,), dtype=numpy.uint64)
n_hist_bins = args.nbins

hist_binbounds = numpy.empty((pcoord_ndim,n_hist_bins+1), dtype=numpy.float64)
hist = numpy.zeros((pcoord_ndim, n_iters, n_hist_bins), numpy.float64)
hist_weights = numpy.zeros((pcoord_ndim, n_iters), numpy.float64)

# Determine range of data
for idim in xrange(0, pcoord_ndim):
    for n_iter in xrange(start_iter, stop_iter+1):
        pcoords = data_manager.get_pcoord_array(n_iter)
        if n_iter == start_iter:
            n_pcoord_points[idim] = pcoords.size
            min_pcoord[idim] = pcoords.min()
            max_pcoord[idim] = pcoords.max()
            #var_pcoord[idim] = pcoords.std()**2        
        else:
            n_pcoord_points[idim] += pcoords.size
            min_pcoord[idim] = min(min_pcoord[idim], pcoords.min())
            max_pcoord[idim] = max(max_pcoord[idim], pcoords.max())
            #var_pcoord[idim] = var_pcoord[idim] + pcoords.std()**2
            
    hist_binbounds[idim] = numpy.linspace(min_pcoord[idim], max_pcoord[idim], n_hist_bins+1)

# Bin data
for idim in xrange(0, pcoord_ndim):
    dx = numpy.diff(hist_binbounds[idim])
    if not args.quiet_mode:
        sys.stdout.write('dimension {}\n'.format(idim))
    for n_iter in xrange(start_iter, stop_iter+1):
        i_iter = n_iter - start_iter
        if not args.quiet_mode and sys.stdout.isatty():
            sys.stdout.write("\r  iteration {}".format(n_iter))
            sys.stdout.flush()
        pcoords = data_manager.get_pcoord_array(n_iter)
        weights = data_manager.get_seg_index(n_iter)['weight']
        
        for iseg in xrange(0, pcoords.shape[0]):
            (lhist, hist_bins) = numpy.histogram(pcoords[iseg], hist_binbounds[idim])
            tweight = weights[iseg] / lhist.sum()
            hist[idim,i_iter,:] += tweight * lhist
            hist_weights[idim,i_iter] += tweight
        
        # Normalize    
        hist[idim,i_iter] /= hist_weights[idim,i_iter]
        I = (hist[idim,i_iter] * dx).sum()
        hist[idim,i_iter] /= I
        assert abs(1 - (hist[idim,i_iter]*dx).sum()) <= 1.0e-15*n_pcoord_points[idim]
        
    if not args.quiet_mode and sys.stdout.isatty():
        sys.stdout.write('\n')
            
    # Plot time course
    if pyplot:
        extent = (min_pcoord[idim], max_pcoord[idim], start_iter, stop_iter)

        if args.timeplot_pattern_lin:
            if args.linear_vrange:
                norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
            else:
                norm = None
            pyplot.figure()
            pyplot.imshow(hist[idim,::-1], aspect='auto', extent=extent, norm=norm, cmap=matplotlib.cm.jet)#, cmap=cm_hovmol)
            pyplot.xlabel('Progress Coordinate (dimension {})'.format(idim))
            pyplot.ylabel('Iteration')
            cb = pyplot.colorbar()
            cb.set_label(r'$P$')
            pyplot.savefig(args.timeplot_pattern_lin % idim)

        if args.timeplot_pattern_log:
            if args.log_vrange:
                lnorm = matplotlib.colors.Normalize(vmin=lvmin, vmax=lvmax)
            else:
                lnorm = None
            lhist = numpy.ma.masked_array(numpy.log10(hist[idim]))
            lhist.mask = ~numpy.isfinite(lhist)
            pyplot.figure()
            pyplot.imshow(lhist[::-1], aspect = 'auto', extent=extent, norm=lnorm, cmap=matplotlib.cm.jet)#, cmap=cm_hovmol)
            pyplot.xlabel('Progress Coordinate (dimension {})'.format(idim))
            pyplot.ylabel('Iteration')
            cb = pyplot.colorbar()
            cb.set_label(r'$\log_{10}\ P$')
            pyplot.savefig(args.timeplot_pattern_log % idim)            
    
    # Now produce average distribution histogram
    # Normalize histogram
    chist = hist[idim].sum(axis=0) / n_iters
    assert abs(1-(chist*dx).sum()) <= 1.0e-15*n_pcoord_points[idim]
    
    # Save histogram
    lb = hist_binbounds[idim, :-1]
    ub = hist_binbounds[idim, 1:]
    midpoints = (ub+lb)/2
    
    if args.output_pattern:
        output_file = open(args.output_pattern % idim, 'wt')
        if not args.suppress_headers:
            output_file.write('# progress coordinate probability distribution for dimension {:d}\n'.format(idim))
            output_file.write('# column 0: lower bound of bin\n')
            output_file.write('# column 1: upper bound of bin\n')
            output_file.write('# column 2: midpoint of bin\n')
            output_file.write('# column 3: probability density\n')
        numpy.savetxt(output_file, numpy.column_stack([lb,ub,midpoints,chist]))
    
    if pyplot and args.plot_pattern:
        pyplot.figure()
        pyplot.plot(midpoints, chist)
        pyplot.xlabel('Progress coordinate (dimension {})'.format(idim))
        pyplot.ylabel('Probability density')
        pyplot.savefig(args.plot_pattern % idim)
        
        
        
