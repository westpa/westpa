# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

import logging
import re, os
from westtools import WESTMasterCommand, WESTSubcommand
import numpy, h5py
from westpa import h5io, textio
import matplotlib
from matplotlib import pyplot
from matplotlib.image import NonUniformImage
from fasthist import normhistnd
from westpa.extloader import get_object

log = logging.getLogger('westtools.plothist')

# Suppress divide-by-zero in log
numpy.seterr(divide='ignore', invalid='ignore')

def sum_except_along(array, axes):
    '''Reduce the given array by addition over all axes except those listed in the scalar or 
    iterable ``axes``'''
    
    try:
        iter(axes)
    except TypeError:
        axes = [axes]
        
    kept = set(axes)
    summed = list(set(range(array.ndim)) - kept)
    
    # Reorder axes so that the kept axes are first, and in the order they 
    # were given
    array = numpy.transpose(array, list(axes) + summed).copy()
    
    # Now, the last len(summed) axes are summed over
    for _ in range(len(summed)):
        array = numpy.add.reduce(array, axis=-1)
        
    return array

class PlotHistBase(WESTSubcommand):
    def __init__(self, parent):
        super(PlotHistBase,self).__init__(parent)
        
        self.input_arg_group = None
        self.output_arg_group = None

        self.input_h5 = None
        self.opmode = None
        self.plotscale = None
        self.enerzero = None
        self.plotrange = None
        self.plottitle = None
        self.postprocess_function = None
        self.plot_contour = None

        # Iteration range for average/evolution
        self.avail_iter_start = None
        self.avail_iter_stop = None
        self.avail_iter_step = None
        self.iter_start = None
        self.iter_stop = None
        self.iter_step = None

        # Iteration for single point
        self.n_iter = None

        # An array of dicts describing what dimensions to work with and
        # what their ranges should be for the plots.
        self.dimensions = []

        self.plot_output_filename = None
        self.text_output_filename = None
        self.hdf5_output_filename = None

    def add_args(self, parser):
        igroup = self.input_arg_group = parser.add_argument_group('input options')
        igroup.add_argument('input', help='HDF5 file containing histogram data')
        igroup.add_argument('firstdim', nargs='?', metavar='DIMENSION',
                            help='''Plot for the given DIMENSION, specified as INT[:[LB,UB]:LABEL], where
                            INT is a zero-based integer identifying the dimension in the histogram, 
                            LB and UB are lower and upper bounds for plotting, and LABEL is the label for
                            the plot axis. (Default: dimension 0, full range.)''')
        
        ogroup = self.output_arg_group = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output', '--plot-output', dest='plot_output', default='hist.pdf', metavar='PLOT_OUTPUT',
                            help='''Store plot as PLOT_OUTPUT. This may be set to an empty string
                            (e.g. --plot-output='') to suppress plotting entirely. The output
                            format is determined by filename extension (and thus defaults to PDF).
                            Default: "%(default)s".''')
        ogroup.add_argument('--hdf5-output',
                            help='''Store plot data in the HDF5 file HDF5_OUTPUT.''')
        ogroup.add_argument('--plot-contour', dest='plot_contour', action='store_const', const=True, default=False,
                            help='''Determines whether or not to superimpose a contour plot over the heatmap for 2D objects.''')
        
        
        pgroup = parser.add_argument_group('plot options')
        pmgroup = pgroup.add_mutually_exclusive_group()
        pgroup.add_argument('--title', dest='title',
                             help='Include TITLE as the top-of-graph title') 
        pmgroup.add_argument('--linear', dest='plotscale', action='store_const', const='linear',
                             help='Plot the histogram on a linear scale.')
        pmgroup.add_argument('--energy', dest='plotscale', action='store_const', const='energy',
                             help='Plot the histogram on an inverted natural log scale, corresponding to (free) energy (default).')
        pmgroup.add_argument('--zero-energy', dest='enerzero', metavar='E', default='min',
                             help='Set the zero of energy to E, which may be a scalar, "min" or "max"')
        pmgroup.add_argument('--log10', dest='plotscale', action='store_const', const='log10',
                             help='Plot the histogram on a base-10 log scale.')
        pgroup.add_argument('--range',
                            help='''Plot histogram ordinates over the given RANGE, specified as "LB,UB",
                            where LB and UB are the lower and upper bounds, respectively. For 1-D plots,
                            this is the Y axis. For 2-D plots, this is the colorbar axis.
                            (Default: full range.)''')
        pgroup.add_argument('--postprocess-function',
                                help='''Names a function (as in module.function) that will be called just prior
                                to saving the plot. The function will be called as ``postprocess(hist, midpoints, binbounds)``
                                where ``hist`` is the histogram that was plotted, ``midpoints`` is the bin midpoints for
                                each dimension, and ``binbounds`` is the bin boundaries for each dimension for 2-D plots,
                                or None otherwise. The plot must be modified in place using the pyplot stateful interface.''')
        
        parser.set_defaults(plotscale='energy')
        
    def process_args(self, args):
        self.plotscale = args.plotscale        
        self.input_h5 = h5py.File(args.input, 'r')
        self.plot_output_filename = args.plot_output
        self.hdf5_output_filename = args.hdf5_output
        self.plot_contour = args.plot_contour
        
        if args.title:
            self.plottitle = args.title
        
        if args.range:
            self.plotrange = self.parse_range(args.range)
        
        if args.firstdim:
            self.dimensions.append(self.parse_dimspec(args.firstdim))
        
        if not args.firstdim:
            self.dimensions.append({'idim': 0, 'label':'dimension 0'})

        if args.enerzero:
            lenerzero = args.enerzero.lower()
            if lenerzero not in ('min', 'max'):
                try:
                    self.enerzero = float(args.enerzero)
                except ValueError:
                    raise ValueError('invalid energy zero point {!r}'.format(args.enerzero))
            else:
                self.enerzero = lenerzero
        else:
            self.enerzero = 'min'

        self.avail_iter_start, self.avail_iter_stop = h5io.get_iter_range(self.input_h5['histograms'])
        try:
            self.avail_iter_step = h5io.get_iter_step(self.input_h5['histograms']) 
        except KeyError:
            self.avail_iter_step = 1
        log.info('HDF5 file {!r} contains data for iterations {} -- {} with a step of {}'.format(args.input,
                                                                                                 self.avail_iter_start,
                                                                                                 self.avail_iter_stop,
                                                                                                 self.avail_iter_step))
        if args.postprocess_function:
            self.postprocess_function = get_object(args.postprocess_function,path=['.'])


    def parse_dimspec(self, dimspec):
        dimdata = {}
        match = re.match(r'([0-9]+)(?::(?:([^,]+),([^:,]+))?(?::(.*))?)?', dimspec)
        if not match:
            raise ValueError('invalid dimension specification {!r}'.format(dimspec))
        
        (idim_txt, lb_txt, ub_txt, label) = match.groups()
        try:
            dimdata['idim'] = int(idim_txt)
            if lb_txt:
                dimdata['lb'] = float(lb_txt)
            if ub_txt:
                dimdata['ub'] = float(ub_txt)
            if label:
                dimdata['label'] = label
            else:
                dimdata['label'] = 'dimension {}'.format(dimdata['idim'])
        except ValueError as e:
            raise ValueError('invalid dimension specification {!r}: {!r}'.format(dimspec, e))
        return dimdata
        
    def parse_range(self, rangespec):
        try:
            (lbt,ubt) = rangespec.split(',')
            return float(lbt), float(ubt)
        except (ValueError,TypeError) as e:
            raise ValueError('invalid range specification {!r}: {!r}'.format(rangespec, e))

    def _ener_zero(self, hist):
        hist = -numpy.log(hist)
        if self.enerzero == 'min':
            hist -= hist.min()
        elif self.enerzero == 'max':
            hist -= hist.max()
        else:
            hist -= self.enerzero
        return hist
        

class PlotSupports2D(PlotHistBase):
    def __init__(self, parent):
        super(PlotSupports2D,self).__init__(parent)
        

    def add_args(self, parser):
        self.input_arg_group.add_argument('seconddim', nargs='?', metavar='ADDTLDIM',
                                          help='''For instantaneous/average plots, plot along the given additional
                                          dimension, producing a color map.''')
        self.output_arg_group.add_argument('--text-output',
                                           help='''Store plot data in a text format at TEXT_OUTPUT. This option is
                                           only valid for 1-D histograms. (Default: no text output.)''')

    def process_args(self, args):
        self.text_output_filename = args.text_output
        if args.seconddim is not None:
            self.dimensions.append(self.parse_dimspec(args.seconddim))

        
    
    def _do_1d_output(self, hist, idim, midpoints):
        enehist = self._ener_zero(hist)
        log10hist = numpy.log10(hist)
        
        if self.hdf5_output_filename:
            with h5py.File(self.hdf5_output_filename, 'w') as output_h5:
                h5io.stamp_creator_data(output_h5)
                output_h5.attrs['source_data'] = os.path.abspath(self.input_h5.filename)
                output_h5.attrs['source_dimension'] = idim
                output_h5['midpoints'] = midpoints
                output_h5['histogram'] = hist
                
        if self.text_output_filename:
            with textio.NumericTextOutputFormatter(self.text_output_filename) as output_file:
                output_file.write_header('source data: {} dimension {}'.format(os.path.abspath(self.input_h5.filename),idim))
                output_file.write_header('column 0: midpoint of bin')
                output_file.write_header('column 1: probability in bin')
                output_file.write_header('column 2: -ln P')
                output_file.write_header('column 3: log10 P')
                numpy.savetxt(output_file, numpy.column_stack([midpoints,hist, enehist, log10hist]))
        
        if self.plot_output_filename:
            if self.plotscale == 'energy':
                plothist = enehist
                label = r'$\Delta F(x)\,/\,kT$' +'\n' + r'$\left[-\ln\,P(x)\right]$'
            elif self.plotscale == 'log10':
                plothist = log10hist
                label = r'$\log_{10}\ P(x)$'
            else:
                plothist = hist
                label = r'$P(x)$'
            pyplot.figure()
            pyplot.plot(midpoints, plothist)
            pyplot.xlim(self.dimensions[0].get('lb'), self.dimensions[0].get('ub'))
            if self.plotrange:
                pyplot.ylim(*self.plotrange)
            pyplot.xlabel(self.dimensions[0]['label'])
            pyplot.ylabel(label)
            if self.plottitle:
                pyplot.title(self.plottitle)
            if self.postprocess_function:
                self.postprocess_function(plothist, midpoints, None)
            pyplot.savefig(self.plot_output_filename)
        
    def _do_2d_output(self, hist, idims, midpoints, binbounds):
        enehist = self._ener_zero(hist)
        log10hist = numpy.log10(hist)
        
        if self.hdf5_output_filename:
            with h5py.File(self.hdf5_output_filename, 'w') as output_h5:
                h5io.stamp_creator_data(output_h5)
                output_h5.attrs['source_data'] = os.path.abspath(self.input_h5.filename)
                output_h5.attrs['source_dimensions'] = numpy.array(idims, numpy.min_scalar_type(max(idims)))
                output_h5.attrs['source_dimension_labels'] = numpy.array([dim['label'] for dim in self.dimensions])
                for idim in idims:
                    output_h5['midpoints_{}'.format(idim)] = midpoints[idim]
                output_h5['histogram'] = hist

                        
        if self.plot_output_filename:
            if self.plotscale == 'energy':
                plothist = enehist
                label = r'$\Delta F(\vec{x})\,/\,kT$' +'\n' + r'$\left[-\ln\,P(x)\right]$'
            elif self.plotscale == 'log10':
                plothist = log10hist
                label = r'$\log_{10}\ P(\vec{x})$'
            else:
                plothist = hist
                plothist[~numpy.isfinite(plothist)] = numpy.nan
                label = r'$P(\vec{x})$'
            
            try:
                vmin, vmax = self.plotrange
            except TypeError:
                vmin, vmax = None, None
                
            pyplot.figure()
            # Transpose input so that axis 0 is displayed as x and axis 1 is displayed as y
#            pyplot.imshow(plothist.T, interpolation='nearest', aspect='auto',
#                          extent=(midpoints[0][0], midpoints[0][-1], midpoints[1][0], midpoints[1][-1]),
#                          origin='lower', vmin=vmin, vmax=vmax)

            # The following reproduces the former calls to imshow and colorbar
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            ax = pyplot.gca()
            nui = NonUniformImage(ax, extent=(midpoints[0][0], midpoints[0][-1], midpoints[1][0], midpoints[1][-1]),
                                  origin='lower', norm=norm)
            nui.set_data(midpoints[0], midpoints[1], plothist.T)
            ax.images.append(nui)
            ax.set_xlim(midpoints[0][0], midpoints[0][-1])
            ax.set_ylim(midpoints[1][0], midpoints[1][-1])
            cb = pyplot.colorbar(nui)
            cb.set_label(label)
            
            pyplot.xlabel(self.dimensions[0]['label'])
            pyplot.xlim(self.dimensions[0].get('lb'), self.dimensions[0].get('ub'))
            pyplot.ylabel(self.dimensions[1]['label'])
            pyplot.ylim(self.dimensions[1].get('lb'), self.dimensions[1].get('ub'))
            if self.plottitle:
                pyplot.title(self.plottitle)
            if self.postprocess_function:
                self.postprocess_function(plothist, midpoints, binbounds)
            if self.plot_contour:
                pyplot.contour(midpoints[0], midpoints[1],plothist.T)
            pyplot.savefig(self.plot_output_filename)


class InstantPlotHist(PlotSupports2D):
    subcommand='instant'
    help_text = 'plot probability distribution for a single WE iteration' 
    description = '''\
Plot a probability distribution for a single WE iteration. The probability
distribution must have been previously extracted with ``w_pdist`` (or, at
least, must be compatible with the output format of ``w_pdist``; see
``w_pdist --help`` for more information).
'''
    def add_args(self, parser):
        self.input_arg_group.add_argument('--iter', metavar='N_ITER', dest='n_iter', type=int,
                                          help='''Plot distribution for iteration N_ITER
                                          (default: last completed iteration).''')        
    def process_args(self, args):
        if args.n_iter:
            self.n_iter = min(args.n_iter, self.avail_iter_stop-1)
        else:
            self.n_iter = self.avail_iter_stop - 1

    def do_instant_plot_1d(self):
        '''Plot the histogram for iteration self.n_iter'''
        
        idim = self.dimensions[0]['idim']
        n_iters = self.input_h5['n_iter'][...]
        iiter = numpy.searchsorted(n_iters, self.n_iter)
        binbounds = self.input_h5['binbounds_{}'.format(idim)][...]
        midpoints = self.input_h5['midpoints_{}'.format(idim)][...]
        hist = self.input_h5['histograms'][iiter]
        
        # Average over other dimensions
        hist = sum_except_along(hist, idim)
        normhistnd(hist, [binbounds])
        self._do_1d_output(hist, idim, midpoints)

    def do_instant_plot_2d(self):
        '''Plot the histogram for iteration self.n_iter'''
        
        idim0 = self.dimensions[0]['idim']
        idim1 = self.dimensions[1]['idim']
        
        n_iters = self.input_h5['n_iter'][...]
        iiter = numpy.searchsorted(n_iters, self.n_iter)
        binbounds_0 = self.input_h5['binbounds_{}'.format(idim0)][...]
        midpoints_0 = self.input_h5['midpoints_{}'.format(idim0)][...]
        binbounds_1 = self.input_h5['binbounds_{}'.format(idim1)][...]
        midpoints_1 = self.input_h5['midpoints_{}'.format(idim1)][...]
        
        hist = self.input_h5['histograms'][iiter]
        
        # Average over other dimensions
        hist = sum_except_along(hist, [idim0,idim1])
        normhistnd(hist, [binbounds_0,binbounds_1])
        self._do_2d_output(hist, [idim0,idim1], [midpoints_0,midpoints_1], [binbounds_0,binbounds_1])

    
    def go(self):
        if len(self.dimensions) == 2:
            self.do_instant_plot_2d()
        else:
            self.do_instant_plot_1d()
        

class AveragePlotHist(PlotSupports2D):
    subcommand='average'
    help_text = 'plot average of a probability distribution over a WE simulation'
    description = '''\
Plot a probability distribution averaged over multiple iterations. The
probability distribution must have been previously extracted with ``w_pdist``
(or, at least, must be compatible with the output format of ``w_pdist``; see
``w_pdist --help`` for more information).
'''
    def add_args(self, parser):
        igroup = self.input_arg_group
        igroup.add_argument('--first-iter', dest='first_iter', type=int, metavar='N_ITER', default=1,
                            help='''Begin averaging at iteration N_ITER (default: %(default)d).''')
        igroup.add_argument('--last-iter', dest='last_iter', type=int, metavar='N_ITER',
                            help='''Conclude averaging with N_ITER, inclusive (default: last completed iteration).''')
        
    def process_args(self, args):
        if args.first_iter:
            self.iter_start = max(args.first_iter, self.avail_iter_start)
        else:
            self.iter_start = self.avail_iter_start
            
        if args.last_iter:
            self.iter_stop = min(args.last_iter+1, self.avail_iter_stop)
        else:
                self.iter_stop = self.avail_iter_stop
        
    
    def do_average_plot_1d(self):
        '''Plot the average histogram for iterations self.iter_start to self.iter_stop'''
        
        idim = self.dimensions[0]['idim']
        n_iters = self.input_h5['n_iter'][...]
        iiter_start = numpy.searchsorted(n_iters, self.iter_start)
        iiter_stop  = numpy.searchsorted(n_iters, self.iter_stop)
        binbounds = self.input_h5['binbounds_{}'.format(idim)][...]
        midpoints = self.input_h5['midpoints_{}'.format(idim)][...]
        #hist = self.input_h5['histograms'][iiter_start:iiter_stop]
        
        for iiter in range(iiter_start, iiter_stop):
            iter_hist = sum_except_along(self.input_h5['histograms'][iiter], idim)
            if iiter == iiter_start:
                hist = iter_hist
            else:
                hist += iter_hist
            del iter_hist

        normhistnd(hist, [binbounds])
        self._do_1d_output(hist, idim, midpoints)


    def do_average_plot_2d(self):
        '''Plot the histogram for iteration self.n_iter'''
        
        idim0 = self.dimensions[0]['idim']
        idim1 = self.dimensions[1]['idim']
        
        n_iters = self.input_h5['n_iter'][...]
        iiter_start = numpy.searchsorted(n_iters, self.iter_start)
        iiter_stop  = numpy.searchsorted(n_iters, self.iter_stop)

        binbounds_0 = self.input_h5['binbounds_{}'.format(idim0)][...]
        midpoints_0 = self.input_h5['midpoints_{}'.format(idim0)][...]
        binbounds_1 = self.input_h5['binbounds_{}'.format(idim1)][...]
        midpoints_1 = self.input_h5['midpoints_{}'.format(idim1)][...]
        
        for iiter in range(iiter_start,iiter_stop): 
            iter_hist = sum_except_along(self.input_h5['histograms'][iiter], [idim0,idim1])
            if iiter == iiter_start:
                hist = iter_hist
            else:
                hist += iter_hist

        normhistnd(hist, [binbounds_0,binbounds_1])
        self._do_2d_output(hist, [idim0,idim1], [midpoints_0, midpoints_1], [binbounds_0,binbounds_1])

    def go(self):
        if len(self.dimensions) == 2:
            self.do_average_plot_2d()
        else:
            self.do_average_plot_1d()
        

class EvolutionPlotHist(PlotHistBase):
    subcommand = 'evolution'
    help_text = 'plot evolution of a probability distribution over the course of a WE simulation'
    description = '''\
Plot a probability distribution as it evolves over iterations. The
probability distribution must have been previously extracted with ``w_pdist``
(or, at least, must be compatible with the output format of ``w_pdist``; see
``w_pdist --help`` for more information).
'''
    def add_args(self, parser):
        igroup = self.input_arg_group
        igroup.add_argument('--first-iter', dest='first_iter', type=int, metavar='N_ITER', default=1,
                            help='''Begin analysis at iteration N_ITER (default: %(default)d).''')
        igroup.add_argument('--last-iter', dest='last_iter', type=int, metavar='N_ITER',
                            help='''Conclude analysis with N_ITER, inclusive (default: last completed iteration).''')
        igroup.add_argument('--step-iter', dest='step_iter', type=int, metavar='STEP',
                            help='''Average in blocks of STEP iterations.''')
        
    def process_args(self, args):
        if args.first_iter:
            self.iter_start = max(args.first_iter, self.avail_iter_start)
        else:
            self.iter_start = self.avail_iter_start
            
        if args.last_iter:
            self.iter_stop = min(args.last_iter+1, self.avail_iter_stop)
        else:
            self.iter_stop = self.avail_iter_stop
            
        if args.step_iter:
            self.iter_step = max(args.step_iter, self.avail_iter_step)
        else:
            self.iter_step = self.avail_iter_step
        log.info('using data for iterations {} -- {} with a step of {}'.format(self.iter_start, self.iter_stop, self.iter_step))
    
    def go(self):
        '''Plot the evolution of the histogram for iterations self.iter_start to self.iter_stop'''
        
        idim = self.dimensions[0]['idim']
        n_iters = self.input_h5['n_iter'][...]
        iiter_start = numpy.searchsorted(n_iters, self.iter_start)
        iiter_stop  = numpy.searchsorted(n_iters, self.iter_stop)
        binbounds = self.input_h5['binbounds_{}'.format(idim)][...]
        midpoints = self.input_h5['midpoints_{}'.format(idim)][...]
        hists_ds = self.input_h5['histograms']
        
        itercount = self.iter_stop - self.iter_start
        
        # We always round down, so that we don't have a dangling partial block at the end
        nblocks = itercount // self.iter_step
            
        block_iters = numpy.empty((nblocks,2), dtype=n_iters.dtype)
        blocked_hists = numpy.zeros((nblocks,hists_ds.shape[1+idim]), dtype=hists_ds.dtype) 
        
        for iblock, istart in enumerate(range(iiter_start, iiter_start+nblocks*self.iter_step, self.iter_step)):
            istop = min(istart+self.iter_step, iiter_stop)
            histslice = hists_ds[istart:istop]
            
            
            # Sum over time
            histslice = numpy.add.reduce(histslice, axis=0)
            
            # Sum over other dimensions
            blocked_hists[iblock] = sum_except_along(histslice, idim)
            
            # Normalize
            normhistnd(blocked_hists[iblock], [binbounds])
            
            block_iters[iblock,0] = n_iters[istart]
            block_iters[iblock,1] = n_iters[istop-1]+1
                        
        #enehists = -numpy.log(blocked_hists)
        enehists = self._ener_zero(blocked_hists)
        log10hists = numpy.log10(blocked_hists)
        
        
        if self.hdf5_output_filename:
            with h5py.File(self.hdf5_output_filename, 'w') as output_h5:
                h5io.stamp_creator_data(output_h5)
                output_h5.attrs['source_data'] = os.path.abspath(self.input_h5.filename)
                output_h5.attrs['source_dimension'] = idim
                output_h5['midpoints'] = midpoints
                output_h5['histograms'] = blocked_hists
                output_h5['n_iter'] = block_iters
                        
        if self.plot_output_filename:
            if self.plotscale == 'energy':
                plothist = enehists
                label = r'$\Delta F(x)\,/\,kT$' +'\n' + r'$\left[-\ln\,P(x)\right]$'
            elif self.plotscale == 'log10':
                plothist = log10hists
                label = r'$\log_{10}\ P(x)$'
            else:
                plothist = blocked_hists
                label = r'$P(x)$'
            
            try:
                vmin, vmax = self.plotrange
            except TypeError:
                vmin, vmax = None, None
                
            pyplot.figure()
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            ax = pyplot.gca()
            nui = NonUniformImage(ax, extent=(midpoints[0], midpoints[-1], block_iters[0,-1], block_iters[-1,-1]),
                                  origin='lower', norm=norm)

            # not sure why plothist works but plothist.T doesn't, and the opposite is true
            # for _do_2d_output
            nui.set_data(midpoints, block_iters[:,-1], plothist)
            ax.images.append(nui)
            ax.set_xlim(midpoints[0], midpoints[-1])
            ax.set_ylim(block_iters[0,-1], block_iters[-1,-1])
            cb = pyplot.colorbar(nui)
            cb.set_label(label)
            pyplot.xlabel(self.dimensions[0]['label'])
            pyplot.xlim(self.dimensions[0].get('lb'), self.dimensions[0].get('ub'))
            pyplot.ylabel('WE Iteration')
            if self.plottitle:
                pyplot.title(self.plottitle)
            if self.postprocess_function:
                self.postprocess_function(plothist, midpoints, binbounds)
            pyplot.savefig(self.plot_output_filename)

    
class PlotHistTool(WESTMasterCommand):
    prog='plothist'
    subparsers_title = 'plotting modes'
    subcommands = [InstantPlotHist,AveragePlotHist,EvolutionPlotHist]
    description = '''\
Plot probability density functions (histograms) generated by w_pdist or other
programs conforming to the same output format. This program operates in one of
three modes:
 
  instant 
    Plot 1-D and 2-D histograms for an individual iteration. See
    ``plothist instant --help`` for more information.
    
  average
    Plot 1-D and 2-D histograms, averaged over several iterations. See
    ``plothist average --help`` for more information.
    
  evolution
    Plot the time evolution 1-D histograms as waterfall (heat map) plots.
    See ``plothist evolution --help`` for more information.
    
This program takes the output of ``w_pdist`` as input (see ``w_pdist --help``
for more information), and can generate any kind of graphical output that
matplotlib supports.


------------------------------------------------------------------------------
Command-line options
------------------------------------------------------------------------------
'''

if __name__ == '__main__':
    PlotHistTool().main()
