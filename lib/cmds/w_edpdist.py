from __future__ import print_function, division; __metaclass__=type
import os, sys, argparse, math, itertools
import numpy
from scipy.stats import kstwobign
import west
from westtools.stats.edfs import EDF

import logging
log = logging.getLogger('w_edpdist')

import westtools
from westtools.aframe import (WESTAnalysisTool,WESTDataReaderMixin,IterRangeMixin,TransitionAnalysisMixin,BinningMixin,
                              KineticsAnalysisMixin,CommonOutputMixin,PlottingMixin)
                              
class WEDPDist(PlottingMixin,CommonOutputMixin,KineticsAnalysisMixin,TransitionAnalysisMixin,
                BinningMixin,IterRangeMixin,WESTDataReaderMixin,WESTAnalysisTool):
    def __init__(self):
        super(WEDPDist,self).__init__()
        self.edp_group = None
        
        self.edfs_group = None
        self.edfs_gname = 'ed_edfs'
        
        self.confidence = None
        self.discard_edf_data = False
        
    def require_edf_group(self):
        self.edfs_group = self.edp_group.require_group(self.edfs_gname)
        for (k,v) in self.get_transitions_ds().attrs.iteritems():
            self.edfs_group.attrs[k] = v
            
    def delete_edf_group(self):
        try:
            del self.edp_group[self.edfs_gname]
        except KeyError:
            pass

    def require_cdfs(self):
        all_pairs = set(self.selected_bin_pair_iter)
        if self.discard_edf_data:
            west.rc.pstatus('Discarding existing EDFs')
            self.delete_edf_group()
            remaining_pairs = all_pairs
        elif self.edfs_gname in self.edp_group:
            avail_pairs = set(self.list_stored_edfs())
            remaining_pairs = all_pairs - avail_pairs

        self.require_edf_group()
        if remaining_pairs:
            west.rc.pstatus('Calculating {:d} EDFs...'.format(len(remaining_pairs)))
            self.calc_cdfs(remaining_pairs)
        else:
            west.rc.pstatus('All required EDFs present.')        
        
    def store_edf_data(self, ibin, fbin, durations, F_all, N_e, ks_probs):
        try:
            del self.edfs_group['{:d}/{:d}'.format(ibin,fbin)]
        except KeyError:
            pass
        
        ibin_group = self.edfs_group.require_group(str(ibin))
        ibin_group.attrs['initial_bin'] = ibin
        fbin_group = ibin_group.require_group(str(fbin))

        fbin_group['durations'] = durations
        
        fbin_group.create_dataset('F_all', data=F_all, compression='gzip', chunks=(4,F_all.shape[1]))
        fbin_group['N_e'] = N_e
        fbin_group['ks_probs'] = ks_probs
        
        for h5obj in (fbin_group, fbin_group['F_all'], fbin_group['N_e'], fbin_group['ks_probs'], fbin_group['durations']):
            h5obj.attrs['dt'] = self.dt
            h5obj.attrs['initial_bin'] = ibin
            h5obj.attrs['final_bin'] = fbin
            self.record_data_iter_range(h5obj)
            self.record_data_iter_step(h5obj)
        
    def list_stored_edfs(self):
        '''Return a list of (ibin,fbin) pairs of EDFs already in the HDF5 file'''
        bin_pairs = []
        for igroup_entry in self.edfs_group:
            try:
                ibin = int(igroup_entry)
            except ValueError:
                continue
            
            for fgroup_entry in self.edfs_group[igroup_entry]:
                try:
                    fbin = int(fgroup_entry)
                except ValueError:
                    continue
                else:
                    bin_pairs.append((ibin,fbin))
        return bin_pairs
    def calc_cdfs(self, bin_pairs):
        transitions_ds = self.get_transitions_ds()
        transitions_niter = transitions_ds['n_iter']
        transitions_ibin = transitions_ds['initial_bin']
        transitions_fbin = transitions_ds['final_bin']
        n_blocks = self.n_iter_blocks()
        last_iters = numpy.empty((n_blocks,), numpy.min_scalar_type(self.last_iter))
        for (ibin,fbin) in bin_pairs:
            tdisc = (transitions_ibin == ibin) & (transitions_fbin == fbin)
            iter_disc = (transitions_niter >= self.first_iter) & (transitions_niter <= self.last_iter)
            transitions = transitions_ds[tdisc & iter_disc]
            if len(transitions):
                # Calculate entire EDF; this must have the highest number of points  
                full_edf = EDF(transitions['duration'] * self.dt, transitions['final_weight'])
                all_durations = full_edf.x
                
                # An array of (n_iter, edf)
                F_all = numpy.empty((n_blocks,len(full_edf)),numpy.float64)
                N_e = numpy.empty((n_blocks,), dtype=numpy.min_scalar_type(len(full_edf)))
                ks_probs = numpy.empty((n_blocks,), numpy.float64)
                
                for iblock, (block_first, block_past_last) in enumerate(self.iter_block_iter()):
                    west.rc.pstatus('\r  {:d}->{:d}: {:d} transitions: iterations [{:d},{:d}]'
                                    .format(ibin,fbin,len(transitions),self.first_iter, block_past_last-1), end='')
                    west.rc.pflush()
                    last_iters[iblock] = block_past_last - 1
                    if iblock == n_blocks - 1:
                        F_all[iblock] = full_edf.F
                        N_e[iblock] = len(full_edf)
                        ks_probs[iblock] = 0
                        break
                     
                    subset_transitions = transitions[transitions['n_iter'] < block_past_last]
                    if len(subset_transitions):
                        subset_edf = EDF(subset_transitions['duration']*self.dt, subset_transitions['final_weight'])
                        N_e[iblock] = len(subset_edf)
                        F_all[iblock] = subset_edf(all_durations)
                        D = numpy.abs(F_all[iblock] - full_edf.F).max()
                        K = math.sqrt(len(F_all)*len(subset_edf)/(len(F_all) + len(subset_edf))) * D
                        if N_e[iblock] < 30:
                            ks_probs[iblock] = 1.0
                        else:
                            ks_probs[iblock] = kstwobign.cdf(K)
                        del subset_edf, subset_transitions
                    else:
                        N_e[iblock] = 0
                        F_all[iblock,:] = 0
                        ks_probs[iblock] = 1

                ks_probs[ks_probs > 1] = 1.0
                ks_probs[ks_probs < 0] = 0.0                       
                self.store_edf_data(ibin, fbin, all_durations, F_all, N_e, ks_probs)
                
                west.rc.pstatus('; N_e = {:d}, mean = {:g}, stdev = {:g}, median = {:g}, 95th percentile = {:g}, max = {:g}'
                                .format(len(full_edf),
                                        float(full_edf.mean()),
                                        float(full_edf.std()),
                                        float(full_edf.quantile(0.5)),
                                        float(full_edf.quantile(0.95)),
                                        float(full_edf.quantile(1.0))), 
                                end='')
            del transitions, all_durations, full_edf
            west.rc.pstatus()

    def write_edfs(self):
        if not self.edf_output_pattern:
            west.rc.pstatus('Not writing any EDFs in text format.')
            return
        else:
            west.rc.pstatus('Writing EDFs to text files...')
        
        bin_labels = list(self.mapper.labels)
        mbw = len(str(len(bin_labels)-1))
        for (ibin, fbin) in self.selected_bin_pair_iter:
            edf_group = self.edfs_group['{:d}/{:d}'.format(ibin,fbin)]
            edf = edf_group['F_all'][-1,:]
            durations = edf_group['durations'][...]
            #ks_probs = edf_group['ks_probs'][...]
            outfile_name = self.edf_output_pattern % (ibin,fbin)
            west.rc.pstatus('  {:{mbw}d}->{:<{mbw}d} to {:s}'.format(ibin,fbin,outfile_name,mbw=mbw))
            outfile = open(outfile_name,'wt')
            
            if not self.output_suppress_headers:
                outfile.write('# transition duration empirical distribution function\n')
                if self.output_print_bin_labels:
                    outfile.write('# initial bin: {:{mbw}d}, {:s}\n'.format(ibin, bin_labels[ibin], mbw=mbw))
                    outfile.write('# final bin:   {:{mbw}d}, {:s}\n'.format(fbin, bin_labels[fbin], mbw=mbw))
                else:
                    outfile.write('# initial bin: {:{mbw}d}\n'.format(ibin, mbw=mbw))
                    outfile.write('# final bin:   {:{mbw}d}\n'.format(fbin, mbw=mbw))
                outfile.write('# transitions collected over iterations [{:d},{:d}]\n'.format(self.first_iter,self.last_iter))
                outfile.write('''\
# --
# column 0: duration
# column 1: P(randomly-chosen duration <= duration from column 0)
''')
    
            numpy.savetxt(outfile, numpy.column_stack([durations, edf]))
            outfile.close()
            del edf, durations, edf_group
    
    def plot_evol(self):
        if not self.evol_plot_pattern:
            west.rc.pstatus('Not plotting evolution of duration distribution.')
            return
        else:
            west.rc.pstatus('Plotting evolution of duration distributions...')
        
        matplotlib = self.require_matplotlib()
        from matplotlib import pyplot
        from matplotlib import gridspec
        
        mbw = len(str(self.n_bins))
        for (ibin,fbin) in self.selected_bin_pair_iter:
            output_filename = self.evol_plot_pattern % (ibin, fbin)
            
            edf_group = self.edfs_group['{:d}/{:d}'.format(ibin,fbin)]
            edfs = self.slice_per_iter_data(edf_group['F_all'])
            durations = edf_group['durations'][...]
            ks_probs = self.slice_per_iter_data(edf_group['ks_probs'])
            
            # Use the smallest difference between measured durations as the resolution for
            # the evolution plot. This can be less than dt due to numerical error, so
            # set to dt if that's the case
            ddiff = numpy.diff(durations)
            dres = max(self.dt, ddiff[ddiff > 0].min())
            interp_durations = numpy.arange(durations[0], durations[-1]+dres, dres)

            # Evaluate the EDF for each iteration at the higher-resolution set of durations;
            # this prepares us to use matplotlib's imshow() to generate the evolution plot
            # while preserving the ED scale            
            interp_edfs = numpy.empty((len(edfs), len(interp_durations)), numpy.float64)
            means = numpy.empty(len(edfs), numpy.float64)
            quantiles = numpy.empty((len(edfs), 5), numpy.float64)
            for iblock in xrange(len(edfs)):
                iter_edf = EDF.from_arrays(durations, edfs[iblock])
                interp_edfs[iblock,:] = iter_edf(interp_durations)
                means[iblock] = iter_edf.mean()
                quantiles[iblock,:] = iter_edf.quantiles([0.05, 0.25, 0.5, 0.75, 0.95])
            
            # Do the plot
            gs = gridspec.GridSpec(2,2, width_ratios = [5,1], height_ratios = [1,15], hspace=0.10, wspace=0.10)
            ax1 = pyplot.subplot(gs[2])
            extent = (durations[0], durations[-1], self.first_iter, self.last_iter)
            im1 = pyplot.imshow(interp_edfs[::-1], aspect='auto', extent=extent, cmap='Spectral')
            
            # poor-man's contour plot, because it's faster
            iter_range = self.iter_range()
            for iquant in xrange(quantiles.shape[-1]):
                #pyplot.plot(quantiles[:,iquant], iter_range, color='white', linewidth=1.5)
                pyplot.plot(quantiles[:,iquant], iter_range, color='black', linewidth=0.75)#, linewidth=0.5)
                
            pyplot.plot(means[:], iter_range, color='black', linestyle='dashed', linewidth=0.75)
                        
            pyplot.xlabel(r'${:d} \rightarrow {:d}$ transition duration $t_\mathrm{{ed}}$'.format(ibin,fbin))
            pyplot.ylabel('Iteration')
            pyplot.xlim(durations[0], durations[-1])
            pyplot.ylim(self.first_iter, self.last_iter)
            
            
            ax2 = pyplot.subplot(gs[3])
            #pyplot.twinx()
            pyplot.plot(numpy.log10(ks_probs), iter_range)
            pyplot.xlim(-6,-4,-2,0)
            pyplot.xticks([-6,0])
            pyplot.vlines(numpy.log10(0.05), iter_range[0], iter_range[-1], linestyles='dashed')
            pyplot.ylim(self.first_iter, self.last_iter)
            pyplot.tick_params(axis='y', labelleft=False, labelright=True)
            pyplot.xlabel(r'$\log_{10}\left(1-P(H_0)\right)$')
            
            ax0 = pyplot.subplot(gs[0])
            cb = pyplot.colorbar(im1, orientation='horizontal', cax=ax0)
            #cb.set_label(r'$F\left(t_\mathrm{ed}\right)$')
            ax0.tick_params(labelbottom=False, labeltop=True, labelsize=8)
            cb.set_ticks([0.0,0.25,0.5,0.75,0.9,0.95,1.0])
            ax0.set_ylabel(r'$F(t_\mathrm{ed})$', rotation=0)
            for label in ax0.get_xticklabels():
                label.set_rotation(90)
                        
            pyplot.savefig(output_filename)
            pyplot.clf()
            west.rc.pstatus('  {:{mbw}d}->{:<{mbw}d} to {:s}'.format(ibin,fbin,output_filename,mbw=mbw))
            
            del ks_probs, durations, interp_durations, interp_edfs, quantiles, means, edfs, edf_group
        
wedp = WEDPDist()

parser = argparse.ArgumentParser('w_edpdist', description='''\
Calculate the transition event duration probability distribution as a function of number of iterations,
and evaluate WEST simulation convergence. 
''')
west.rc.add_args(parser)
wedp.add_args(parser)

edfgroup = parser.add_argument_group('EDF options')
edfgroup.add_argument('--discard-edf-data', action='store_true',
                      help='''Discard any existing EDF data for event durations''')

ogroup = parser.add_argument_group('output options')
ogroup.add_argument('--summary-output-pattern', default='edf_ed_summary_%d_%d.txt',
                    help='''Write time-ordered summary statistics (mean, median, 95th percentile, etc.) to files named according to
                    SUMMARY_OUTPUT_PATTERN, which must contain two printf-style escape sequences which will be replaced
                    with the indices of the initial and final bins whose duration EDF is being considered
                    (default: %(default)s).''')
ogroup.add_argument('--edf-output-pattern', default='edf_ed_%d_%d.txt',
                    help='''Write entire (final) duration EDF to files named according to EDF_OUTPUT_PATTERN, which must
                    contain two printf-style escape sequences which will be replaced with the indices of the initial and
                    final bins whose duration EDF is being considered (default: %(default)s).''')

pgroup = parser.add_argument_group('plotting options')
pgroup.add_argument('--evol-plot-pattern', default='edf_ed_evol_%d_%d.pdf',
                    help='''Plot EDF as a function of simulation time to files named according to EVOL_PLOT_PATTERN,
                    which must contain two printf-style escape sequences which will be replaced with the indices of the
                    initial and final bins whose duration EDF is being considered (default: %(default)s).''')

wedp.add_common_output_args(ogroup)

args = parser.parse_args()

west.rc.process_args(args, config_required=False)
wedp.process_args(args)
wedp.process_common_output_args(args)
wedp.dt = args.dt
wedp.discard_edf_data = args.discard_edf_data
wedp.summary_output_pattern = args.summary_output_pattern
wedp.edf_output_pattern = args.edf_output_pattern
wedp.evol_plot_pattern = args.evol_plot_pattern
for pattern in (args.summary_output_pattern, args.edf_output_pattern, args.evol_plot_pattern):
    if pattern:
        try:
            pattern % (0,0)
        except TypeError:
            raise ValueError('invalid pattern {!r}'.format(pattern))
wedp.suppress_headers = args.suppress_headers
wedp.print_bin_labels = args.print_bin_labels

# Preliminary checks and required calculation steps
wedp.check_iter_range()
wedp.check_bin_selection()
wedp.open_analysis_backing()
wedp.edp_group = wedp.require_analysis_group('w_edpdist')
wedp.require_bin_assignments()
wedp.require_transitions()
wedp.require_edf_group()
wedp.require_cdfs()
wedp.write_edfs()
if wedp.matplotlib_avail:
    wedp.plot_evol()
else:
    west.rc.pstatus('matplotlib not available; not generating plots')

