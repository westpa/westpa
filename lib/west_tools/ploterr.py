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

from __future__ import print_function, division; __metaclass__ = type
import logging
import re, os
from westtools import WESTMasterCommand, WESTSubcommand, ProgressIndicatorComponent
import numpy, h5py
from westpa import h5io
import matplotlib
from matplotlib import pyplot

log = logging.getLogger('westtools.ploterrs')

class CommonPloterrs(WESTSubcommand):
    def __init__(self, parent):
        super(CommonPloterrs,self).__init__(parent)
        
        self.progress = ProgressIndicatorComponent()

        self.xscale = None
        self.yscale = None
        self.xrange = None
        self.yrange = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        
        self.plot_options_group = None
        
    def add_args(self, parser):
        self.progress.add_args(parser)

        pogroup = self.plot_options_group = parser.add_argument_group('plot options')
        pogroup.add_argument('--xscale', choices=['linear', 'log', 'symlog'], default='linear',
                             help='''Use "linear", "log", or "symlog" scaling for the x axis.
                             (Default: %(default)s).''')
        pogroup.add_argument('--yscale', choices=['linear', 'log', 'symlog'], default='linear',
                             help='''Use "linear", "log", or "symlog" scaling for the y axis.
                             (Default: %(default)s).''')
        pogroup.add_argument('--xrange',
                             help='''Restrict X range to XRANGE, which must be formatted as "xmin,xmax".
                             (Default: determined by input data.)''')
        pogroup.add_argument('--yrange',
                             help='''Restrict Y range to YRANGE, which must be formatted as "ymin,ymax".
                             (Default: determined by input data.)''')
        pogroup.add_argument('--xlabel',
                             help='''Use XLABEL for the x-axis label. (Default: varies.)''')
        pogroup.add_argument('--ylabel',
                             help='''Use YLABEL for the y-axis label. (Default: varies.)''')
        pogroup.add_argument('--title',
                             help='''Use TITLE for the plot title. (Default: varies.)''')

    def process_args(self, args):
        self.progress.process_args(args)
        
        if args.xrange:
            self.xrange = self.parse_range(args.xrange)

        if args.yrange:
            self.yrange = self.parse_range(args.yrange)
            
        self.xscale = args.xscale
        self.yscale = args.yscale
        self.xlabel = args.xlabel or 'Iteration'
        self.ylabel = args.ylabel
        self.title = args.title

    def parse_range(self, rangespec):
        try:
            (lbt,ubt) = rangespec.split(',')
            return float(lbt), float(ubt)
        except (ValueError,TypeError) as e:
            raise ValueError('invalid range specification {!r}: {!s}'.format(rangespec, e))

    def do_plot(self, data, output_filename, title=None, x_range=None, y_range=None, x_label=None, y_label=None):
        if not output_filename:
            return
        
        title = title or self.title
        x_range = x_range or self.xrange
        y_range = y_range or self.yrange
        x_label = x_label or self.xlabel
        y_label = y_label or self.ylabel 
        
        iters = data['iter_stop'] - 1

        pyplot.figure()
        pyplot.plot(iters, data['expected'], color='black')
        pyplot.plot(iters, data['ci_lbound'], color='gray')
        pyplot.plot(iters, data['ci_ubound'], color='gray')
        
        pyplot.gca().set_xscale(self.xscale)
        pyplot.gca().set_yscale(self.yscale)
        
        if title:
            pyplot.title(title)
        
        if x_range is not None:
            pyplot.xlim(x_range)
        
        if y_range is not None:
            pyplot.ylim(y_range)
        
        if x_label:
            pyplot.xlabel(x_label)
        
        if y_label:
            pyplot.ylabel(y_label)

        pyplot.savefig(output_filename)


class GenericIntervalSubcommand(CommonPloterrs):
    description = '''\
Plots generic expectation/CI data. A path to the HDF5 file and the dataset
within it must be provided. This path takes the form **FILENAME/PATH[SLICE]**.
If the dataset is not a vector (one dimensional) then a slice must be provided.
For example, to access the state 0 to state 1 rate evolution calculated by
``w_kinavg``, one would use ``kinavg.h5/rate_evolution[:,0,1]``.


-----------------------------------------------------------------------------
Command-line arguments
-----------------------------------------------------------------------------
'''
    subcommand = 'generic'
    help_text = 'arbitrary HDF5 file and dataset'
    def __init__(self, parent):
        super(GenericIntervalSubcommand,self).__init__(parent)
        self.h5file = None
        self.h5dset = None
        self.dset_slice = None
        self.output_filename = None

    def add_args(self, parser):
        iogroup = parser.add_argument_group('input/output options')
        iogroup.add_argument('-o', '--output', default='errbars.pdf',
                             help='''Write plot to OUTPUT (default: %(default)s), whose format will
                             be determined by filename extension.''')
        iogroup.add_argument('dsspec',
                             help='''Use data located at DSSPEC, which must be formatted as 
                             FILENAME/PATH[SLICE]. FILENAME is the HDF5 file to read, PATH is the
                             HDF5 path to the dataset, and SLICE, if provided, must be the Numpy-style
                             slice (including brackets) which selects a vector of data of the
                             appropriate type.''')
        
    def process_args(self, args):
        self.output_filename = args.output
        (pathname, slicestr) = re.search(r'([^[]+)(\[[^\]]+\])?$', args.dsspec).groups()
        if slicestr:
            sl = eval('numpy.index_exp' + slicestr)
        else:
            sl = numpy.index_exp[...]
        self.h5file, self.h5dset = h5io.resolve_filepath(pathname, mode='r')
        self.dset_slice = sl
        
        
    def load_and_validate_data(self):
        reqd_fields = set(['iter_start', 'iter_stop', 'expected', 'ci_lbound', 'ci_ubound'])
        
        self.progress.indicator.new_operation('loading data')
        data = self.h5dset[self.dset_slice]
        
        if data.ndim != 1:
            raise TypeError('dataset to be plotted must be 1-dimensional')
        try:
            fieldnames = set(data.dtype.fields.keys())
        except AttributeError:
            raise TypeError('dataset has inappropriate type')
        else:
            if len(fieldnames & reqd_fields) < len(reqd_fields):
                raise TypeError('dataset does not contain correct fields')
            
        return data
    

    def go(self):
        with self.progress.indicator:
            data = self.load_and_validate_data()
            self.progress.indicator.new_operation('plotting')
            self.do_plot(data, self.output_filename)

class DirectKinetics(CommonPloterrs):
    subcommand = 'd.kinetics'
    help_text = 'output of w_direct kinetics'
    input_filename = 'direct.h5'
    flux_output_filename = 'flux_evolution_d_{state_label}.pdf'
    rate_output_filename = 'rate_evolution_d_{istate_label}_{fstate_label}.pdf'
    description = '''\
Plot evolution of state-to-state rates and total flux into states as generated
by ``w_{direct/reweight} kinetics`` (when used with the ``--evolution-mode`` 
option). Plots are generated for all rates/fluxes calculated. Output filenames 
require (and plot titles and axis labels support) substitution based on which 
flux/rate is being plotted:

  istate_label, fstate_label
    *(String, for rates)* Names of the initial and final states, as originally
    given to ``w_assign``.

  istate_index, fstate_index
    *(Integer, for rates)* Indices of initial and final states.
    
  state_label
    *(String, for fluxes)* Name of state
    
  state_index
    *(Integer, for fluxes)* Index of state
'''
    
    def __init__(self, parent):
        super(DirectKinetics,self).__init__(parent)
        self.kinavg_file = None

        self.dset_slice = None
        self.rate_output_pattern = None
        self.flux_output_pattern = None
        
        self.state_labels = None


    def add_args(self, parser):
        iogroup = parser.add_argument_group('input/output')
        iogroup.add_argument('-i', '--input', default=self.input_filename,
                             help='''Read kinetics results from INPUT (default: %(default)s).''')
        iogroup.add_argument('--rate-output', default=self.rate_output_filename,
                             help='''Filename pattern for rate evolution output. See above for valid
                             field names. (Default: %(default)r).''')
        iogroup.add_argument('--flux-output', default=self.flux_output_filename,
                             help='''Filename pattern for flux evolution output. See above for valid
                             field names. (Default: %(default)r).''')
    
    def process_args(self, args):
        self.kinavg_file = h5py.File(args.input, 'r')
        
        self.state_labels = list(self.kinavg_file['state_labels'][...])
        
        self.rate_output_pattern = args.rate_output
        self.flux_output_pattern = args.flux_output
                
    def plot_flux(self, istate):
        label = self.state_labels[istate]
        data = self.kinavg_file['target_flux_evolution'][:,istate]
        
        if (data['iter_start'] == 0).all():
            # No data
            return
        
        subdict = dict(state_label=label, state_index=istate)
        
        output_filename = self.flux_output_pattern.format(**subdict) if self.flux_output_pattern else None
        
        title = self.title if self.title is not None else 'Flux into state "{state_label}"'
        title = title.format(**subdict)
        
        x_label = self.xlabel.format(**subdict) if self.xlabel else None
        
        y_label = self.ylabel if self.ylabel is not None else r'Flux $(\tau^{{-1}})$'
        y_label = y_label.format(**subdict)
        
        self.do_plot(data, output_filename, title, x_label=x_label, y_label=y_label)
    
    def plot_rate(self, istate, jstate):
        ilabel = self.state_labels[istate]
        jlabel = self.state_labels[jstate]
        data = self.kinavg_file['rate_evolution'][:,istate, jstate]
        
        if (data['iter_start'] == 0).all():
            # No data
            return
        
        subdict = dict(istate_label=ilabel, istate_index=istate,
                       fstate_label=jlabel, fstate_index=jstate)
        
        output_filename = self.rate_output_pattern.format(**subdict) if self.rate_output_pattern else None

        title = self.title if self.title is not None else 'Rate from state "{istate_label}" to state "{fstate_label}"'
        title = title.format(**subdict)
        
        x_label = self.xlabel.format(**subdict) if self.xlabel else None
        
        y_label = self.ylabel if self.ylabel is not None else r'Rate $(\tau^{{-1}})$'
        y_label = y_label.format(**subdict)
        
        self.do_plot(data, output_filename, title, x_label=x_label, y_label=y_label)
    
    def go(self):
        pi = self.progress.indicator
        nstates = len(self.state_labels)
        with pi:
            # if --evolution-mode wasn't specified, neither of these exist:
            if 'target_flux_evolution' in self.kinavg_file:
                pi.new_operation('plotting fluxes', nstates)
                for istate in xrange(nstates):
                    self.plot_flux(istate)
                    pi.progress += 1
            
            # if --evolution-mode wasn't specified, we won't get this either
            if 'rate_evolution' in self.kinavg_file:
                pi.new_operation('plotting rates', nstates*nstates)
                for istate in xrange(nstates):
                    for jstate in xrange(nstates):
                        self.plot_rate(istate, jstate)
                        pi.progress += 1
            else:
                print('rate evolution not available')

class DirectStateprobs(CommonPloterrs):
    subcommand = 'd.probs'
    help_text = 'output of w_direct probs'
    input_filename = 'direct.h5'
    pop_output_filename = 'pop_evolution_d_{state_label}.pdf'
    color_output_filename = 'color_evolution_d_{state_label}.pdf'
    description = '''\
Plot evolution of macrostate populations and associated uncertainties. Plots
are generated for all states calculated. Output filenames require (and plot
titles and axis labels support) substitution based on which state is being
plotted:

  state_label
    *(String, for fluxes)* Name of state
    
  state_index
    *(Integer, for fluxes)* Index of state
'''
    
    def __init__(self, parent):
        super(DirectStateprobs,self).__init__(parent)
        self.stateprobs_file = None

        self.dset_slice = None
        self.rate_output_pattern = None
        self.flux_output_pattern = None
        
        self.state_labels = None


    def add_args(self, parser):
        iogroup = parser.add_argument_group('input/output')
        iogroup.add_argument('-i', '--input', default=self.input_filename,
                             help='''Read w_kinavg results from INPUT (default: %(default)s).''')
        iogroup.add_argument('--population-output', default=self.pop_output_filename,
                             help='''Filename pattern for population evolution output. See above for valid
                             field names. (Default: %(default)r).''')
        iogroup.add_argument('--color-output', default=self.color_output_filename,
                             help='''Filename pattern for ensemble evolution output. See above for valid
                             field names. (Default: %(default)r).''')
    
    def process_args(self, args):
        self.stateprobs_file = h5py.File(args.input, 'r')
        
        self.state_labels = list(self.stateprobs_file['state_labels'][...])
        self.pop_output_pattern = args.population_output
        self.color_output_pattern = args.color_output
                
    def plot_pop(self, istate):
        label = self.state_labels[istate]
        data = self.stateprobs_file['state_pop_evolution'][:,istate]
        
        if (data['iter_start'] == 0).all():
            # No data
            return
        
        subdict = dict(state_label=label, state_index=istate)
        
        output_filename = self.pop_output_pattern.format(**subdict) if self.pop_output_pattern else None
        
        title = self.title if self.title is not None else 'Population in state "{state_label}"'
        title = title.format(**subdict)
        
        x_label = self.xlabel.format(**subdict) if self.xlabel else None
        y_label = self.ylabel if self.ylabel is not None else r'Population'
        y_label = y_label.format(**subdict)
        
        self.do_plot(data, output_filename, title, x_label=x_label, y_label=y_label)

    def plot_color(self, istate):
        label = self.state_labels[istate]
        data = self.stateprobs_file['color_prob_evolution'][:,istate]
        
        if (data['iter_start'] == 0).all():
            # No data
            return
        
        subdict = dict(state_label=label, state_index=istate)
        
        output_filename = self.color_output_pattern.format(**subdict) if self.color_output_pattern else None
        
        title = self.title if self.title is not None else 'Population in ensemble "{state_label}"'
        title = title.format(**subdict)
        
        x_label = self.xlabel.format(**subdict) if self.xlabel else None
        y_label = self.ylabel if self.ylabel is not None else r'Population'
        y_label = y_label.format(**subdict)
        
        self.do_plot(data, output_filename, title, x_label=x_label, y_label=y_label)
    
    def go(self):
        pi = self.progress.indicator
        nstates = len(self.state_labels)
        with pi:
            if 'state_pop_evolution' in self.stateprobs_file:
                pi.new_operation('plotting populations', nstates)
                for istate in xrange(nstates):
                    self.plot_pop(istate)
                    pi.progress += 1

            if 'color_prob_evolution' in self.stateprobs_file:
                pi.new_operation('plotting ensemble populations', nstates)
                for istate in xrange(nstates):
                    self.plot_color(istate)
                    pi.progress += 1
            else:
                print('population evolution not available')

class ReweightStateprobs(DirectStateprobs):
    subcommand = 'rw.probs'
    help_text = 'output of w_reweight probs'
    input_filename = 'reweight.h5'
    pop_output_filename = 'pop_evolution_rw_{state_label}.pdf'
    color_output_filename = 'color_evolution_rw_{state_label}.pdf'

class ReweightKinetics(DirectKinetics):
    subcommand = 'rw.kinetics'
    help_text = 'output of w_reweight kinetics'
    input_filename = 'reweight.h5'
    flux_output_filename = 'flux_evolution_rw_{state_label}.pdf'
    rate_output_filename = 'rate_evolution_rw_{istate_label}_{fstate_label}.pdf'
            
class PloterrsTool(WESTMasterCommand):
    prog='ploterrs'
    subcommands = [DirectKinetics,DirectStateprobs,ReweightStateprobs,ReweightKinetics,GenericIntervalSubcommand]
    subparsers_title = 'supported input formats'
    description = '''\
Plots error ranges for weighted ensemble datasets. 


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------
'''

if __name__ == '__main__':
    PloterrsTool().main()



