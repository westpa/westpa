from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, itertools, warnings
import numpy
import west

import logging
log = logging.getLogger('w_pcpdist')

import westtools
from westtools.aframe import (WESTAnalysisTool,WESTDataReaderMixin,IterRangeMixin,CommonOutputMixin,PlottingMixin)
                              
class WPCPDist(PlottingMixin,CommonOutputMixin,
                 IterRangeMixin,WESTDataReaderMixin,WESTAnalysisTool):

    def __init__(self):
        super(WPCPDist,self).__init__()
        
        self.pcp_group = None
        
        self.nbins = None
        self.plot_pattern = None
        self.timeplot_pattern_lin = None
        self.timeplot_pattern_log = None
        self.vmin = None
        self.vmax = None
        self.lvmin = None
        self.lvmax = None
        
        self.pcoord_dtype = None
        self.pcoord_ndim = None
        
        self.dims = None
        
        self.include_args['IterRangeMixin']['iter_step'] = False
        
        
    def require_pdists(self):
        self.pcoord_dtype = self.get_pcoord_dataset(self.first_iter).dtype
        self.pcoord_ndim = self.get_pcoord_dataset(self.first_iter).shape[2]
        
        if not self.dims:
            self.dims = list(xrange(self.pcoord_ndim))
        
        for idim in self.dims:
            # Does the group exist?
            # If so, does it bear data for the correct set of iterations?
            # If not, calculate
            groupname = 'dim_{:02d}'.format(idim)
            try:
                is_valid = self.check_data_iter_range_least(self.pcp_group[groupname])
            except KeyError:
                is_valid = False
            if not is_valid:
                try:
                    del self.pcp_group[groupname]
                except KeyError:
                    pass
                
                dimgroup = self.pcp_group.create_group(groupname)
                self.record_data_iter_range(dimgroup)
                self.calc_pdist(idim, dimgroup)
            else:
                west.rc.pstatus('histogram for dimension {:d} already prepared'.format(idim))
                
    def do_plots(self):
        for idim in self.dims:
            groupname = 'dim_{:02d}'.format(idim)
            dimgroup = self.pcp_group[groupname]
            if self.timeplot_pattern_lin:
                ofn = self.timeplot_pattern_lin % idim
                self.plot_evol_lin(idim, dimgroup, ofn)
            if self.timeplot_pattern_log:
                ofn = self.timeplot_pattern_log % idim
                self.plot_evol_log(idim, dimgroup, ofn)
            if self.plot_pattern:
                ofn = self.plot_pattern % idim
                self.plot_pdist(idim, dimgroup, ofn)
                    
    def calc_pdist(self, idim, dimgroup):
        west.rc.pstatus('calculating histograms for dimension {:d}'.format(idim))
        
        min_pcoord = None
        max_pcoord = None
        n_points = None
        n_iters = self.last_iter - self.first_iter + 1
        
        hist = numpy.zeros((n_iters,self.nbins), dtype=numpy.float64)
        hist_weights = numpy.zeros((n_iters,), dtype=numpy.float64)
        
        # Determine range of data
        for n_iter in xrange(self.first_iter,self.last_iter+1):
            pcoord_ds = self.get_pcoord_dataset(n_iter)
            pcoords = pcoord_ds[:,:,idim]
            if min_pcoord is None:
                min_pcoord = pcoords.min()
                max_pcoord = pcoords.max()
                n_points = len(pcoords)
            else:
                min_pcoord = min(min_pcoord, pcoords.min())
                max_pcoord = max(max_pcoord, pcoords.max())
                n_points += len(pcoords)
            del pcoords, pcoord_ds
                
        hist_binbounds = numpy.linspace(min_pcoord, max_pcoord, self.nbins+1)
        dimgroup['binbounds'] = hist_binbounds
        
        # Bin data
        dx = numpy.diff(hist_binbounds)
        for n_iter in xrange(self.first_iter,self.last_iter+1):
            iiter = n_iter - self.first_iter
            pcoords = self.get_pcoord_dataset(n_iter)[:,:,idim]
            weights = self.get_seg_index(n_iter)['weight']
            for iseg in xrange(pcoords.shape[0]):
                (lhist, hist_bins) = numpy.histogram(pcoords[iseg], hist_binbounds)
                tweight = weights[iseg] / lhist.sum()
                hist[iiter] += tweight * lhist
                hist_weights[iiter] += tweight
                
                del tweight, lhist, hist_bins
                
            # Normalize
            hist[iiter] /= hist_weights[iiter]
            I = (hist[iiter]*dx).sum()
            hist[iiter] /= I
            
            del pcoords, weights
            
        dimgroup['histogram'] = hist
        self.record_data_iter_range(dimgroup['histogram'])
             
    def plot_evol_lin(self, idim, dimgroup, ofn):
        west.rc.pstatus('Saving time evolution of dimension {:d} density (linear scale) to "{:s}"'.format(idim, ofn))
        matplotlib = self.require_matplotlib()
        from matplotlib import pyplot
        binbounds = dimgroup['binbounds'][:]
        hist=self.slice_per_iter_data(dimgroup['histogram'])
        extent = (binbounds[0], binbounds[-1], self.first_iter, self.last_iter)
        if self.vmin is not None:
            norm = matplotlib.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
        else:
            norm = None
            
        pyplot.figure()
        pyplot.imshow(hist[::-1], aspect='auto', extent=extent, norm=norm, cmap=matplotlib.cm.jet)
        pyplot.xlabel('Progress Coordinate (dimension {:d})'.format(idim))
        pyplot.ylabel('Iteration')
        cb = pyplot.colorbar()
        cb.set_label(r'$P$')
        pyplot.savefig(ofn)
        
    def plot_evol_log(self, idim, dimgroup, ofn):
        west.rc.pstatus('Saving time evolution of dimension {:d} density (log scale) to "{:s}"'.format(idim, ofn))
        matplotlib = self.require_matplotlib()
        from matplotlib import pyplot
        binbounds = dimgroup['binbounds'][:]
        hist=self.slice_per_iter_data(dimgroup['histogram'])
        extent = (binbounds[0], binbounds[-1], self.first_iter, self.last_iter)
        if self.lvmin is not None:
            norm = matplotlib.colors.Normalize(vmin=self.lvmin, vmax=self.lvmax)
        else:
            norm = None
            
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            lhist = numpy.ma.masked_array(numpy.log10(hist))
        lhist.mask = ~numpy.isfinite(lhist)
            
        pyplot.figure()
        pyplot.imshow(lhist[::-1], aspect='auto', extent=extent, norm=norm, cmap=matplotlib.cm.jet)
        pyplot.xlabel('Progress Coordinate (dimension {:d})'.format(idim))
        pyplot.ylabel('Iteration')
        cb = pyplot.colorbar()
        cb.set_label(r'$\log_{10}\ P$')
        pyplot.savefig(ofn)
        
    def plot_pdist(self, idim, dimgroup, ofn):
        west.rc.pstatus('Plotting average distribution of dimension {:d} to "{:s}"'.format(idim, ofn))
        matplotlib = self.require_matplotlib()
        from matplotlib import pyplot
        binbounds = dimgroup['binbounds'][:]
        hist=self.slice_per_iter_data(dimgroup['histogram'])
        ahist = numpy.mean(hist, axis=0)
        midpoints = (binbounds[1:] + binbounds[:-1]) / 2
        pyplot.figure()
        pyplot.plot(midpoints, ahist)
        pyplot.xlabel('Progress coordinate (dimension {:d})'.format(idim))
        pyplot.ylabel('Probability density')
        pyplot.savefig(ofn)
        
wpcp = WPCPDist()

parser = argparse.ArgumentParser('w_pcpdist', description='''\
Calculate and plot probability distributions of progress coordinate values. 
''')
west.rc.add_args(parser)
wpcp.add_args(parser)

hgroup = parser.add_argument_group('histogram options')
hgroup.add_argument('-B', '--nbins', dest='nbins', type=int, default=1000,
                    help='''Use NBINS for each progress coordinate histogram (default: %(default)s).''')
hgroup.add_argument('--discard-pcpdata', action='store_true',
                    help='''Discard any existing progress coordinate probability distribution data.''')
hgroup.add_argument('dims', metavar='DIM', nargs='*', type=int,
                    help='''Calculate progress coordinate probability distribution(s) for the given dimension(s)
                    (default: calculate for all dimensions)''')

ogroup = parser.add_argument_group('output options')
ogroup.add_argument('--plot', dest='plot_pattern', default='pcprob_%d.pdf',
                    help='Store a plot of the probability distribution of each progress coordinate dimension '
                        +'in PLOT_PATTERN, which must contain a formatting code which will be '
                        +'replaced by the progress coordinate dimension (default: pcprob_%%d.pdf)')
ogroup.add_argument('--timeplot-lin', dest='timeplot_pattern_lin', default='pcprob_evol_%d_lin.pdf',
                    help='Store a colormap plot of the evolution of the probability distribution of '
                        +'each coordinate dimension in TIMEPLOT_PATTERN_LIN, which must contain a formatting ' 
                        +'code which will be replaced by the progress coordinate dimension  '
                        +'(default: pcprob_evol_%%d_lin.pdf)')
ogroup.add_argument('--timeplot-log', dest='timeplot_pattern_log', default='pcprob_evol_%d_log.pdf',
                    help='Store a colormap plot of the evolution of the log of the probability distribution of '
                        +'each coordinate dimension in TIMEPLOT_PATTERN_LOG, which must contain a formatting ' 
                        +'code which will be replaced by the progress coordinate dimension  '
                        +'(default: pcprob_evol_%%d_log.pdf)')    
ogroup.add_argument('--vrange', dest='linear_vrange',
                    help='Use LINEAR_VRANGE (a comma-separated pair of values) as the range for the color bar '
                        +'for the linear-scale probability evolution plot. (Default: range of probability density.)')
ogroup.add_argument('--lvrange', dest='log_vrange',
                    help='Use LOG_VRANGE (a comma-separated pair of values) as the range for the color bar '
                        +'for the log-scale probability evolution plot. (Default: range of log of nonzero probability density)')    
wpcp.add_common_output_args(ogroup)

args = parser.parse_args()
west.rc.process_args(args, config_required=False)
wpcp.process_args(args)
wpcp.process_common_output_args(args)

wpcp.nbins = args.nbins
wpcp.plot_pattern = args.plot_pattern
wpcp.timeplot_pattern_lin = args.timeplot_pattern_lin
wpcp.timeplot_pattern_log = args.timeplot_pattern_log
try:
    wpcp.vmin, wpcp.vmax = map(float,args.linear_vrange.split(','))
except (ValueError, IndexError, TypeError):
    sys.stderr.write('invalid argument for --vrange: {!r}\n'.format(args.linear_vrange))
except AttributeError:
    wpcp.vmin = None
    wpcp.vmax = None
try:
    wpcp.lvmin, wpcp.lvmax = map(float,args.log_vrange.split(','))
except (ValueError, IndexError, TypeError):
    sys.stderr.write('invalid argument for --lvrange: {!r}\n'.format(args.log_vrange))
except AttributeError:
    wpcp.lvmin = None
    wpcp.lvmax = None
wpcp.check_iter_range()
wpcp.open_analysis_backing()
if args.discard_pcpdata:
    try:
        del wpcp.anal_h5file['w_pcpdist']
    except KeyError:
        pass
wpcp.pcp_group = wpcp.require_analysis_group('w_pcpdist')
wpcp.dims = args.dims
wpcp.require_pdists()
wpcp.do_plots()
