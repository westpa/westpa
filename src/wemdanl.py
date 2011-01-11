from __future__ import division

import os, sys
from optparse import OptionParser
import wemd
from wemd import Segment, WESimIter

from wemd.util.wetool import WECmdLineMultiTool
from wemd.rc import EX_ERROR, EX_USAGE_ERROR

from logging import getLogger
log = getLogger(__name__)
    
class WEMDAnlTool(WECmdLineMultiTool):
    def __init__(self):
        super(WEMDAnlTool,self).__init__()
        cop = self.command_parser
        
        cop.add_command('lstrajs',
                        'list trajectories',
                        self.cmd_lstrajs, True)
        cop.add_command('transanl',
                        'analyze transitions',
                        self.cmd_transanl, True)
        cop.add_command('fluxanl',
                        'analyze probability flux',
                        self.cmd_fluxanl, True)
        cop.add_command('tracetraj',
                        'trace an individual trajectory and report',
                        self.cmd_tracetraj, True)
        self.add_rc_option(cop)
        
    def add_iter_param(self, parser):
        parser.add_option('-i', '--iteration', dest='we_iter', type='int',
                          help = 'use trajectories as of iteration number '
                               + 'WE_ITER (default: the last complete '
                               + ' iteration)'
                         )
        
    def get_sim_iter(self, we_iter, complete = False):
        self.load_sim_manager()
        self.sim_manager.load_data_manager()
                
        data_manager = self.sim_manager.data_manager
        if we_iter is not None:
            self.we_iter = data_manager.get_we_sim_iter( we_iter )
        else:
            if complete:
                self.we_iter = data_manager.get_last_complete_iter()
            else:
                self.we_iter = data_manager.get_last_iter()
        
        return self.we_iter
                
    def cmd_lstrajs(self, args):
        parser = self.make_parser(description = 'list trajectories')
        self.add_iter_param(parser)
        parser.add_option('-t', '--type', dest='traj_type',
                          type='choice', choices=('live', 'complete', 
                                                  'merged', 'recycled', 'all'),
                          default='live',
                          help='''list all trajectories ("all"), those that are 
                               alive ("live", default), those that are complete
                               for any reason ("complete"), those that have been
                               terminated because of a merge ("merged"),
                               or those that have been terminated and
                               recycled ("recycled")''')
        (opts,args) = parser.parse_args(args)
        
        self.load_sim_manager()
        self.sim_manager.load_data_manager()
        
        data_manager = self.sim_manager.data_manager
        we_iter = self.get_sim_iter(opts.we_iter)
        from wemd.core import Segment
        
        sys.stdout.write('# %s trajectories as of iteration %d:\n'
                                 % (opts.traj_type, we_iter.n_iter))
        sys.stdout.write('#%-11s    %-12s    %-21s\n'
                                 % ('seg_id', 'n_iter', 'weight'))
        
        if opts.traj_type == 'live':
            segments = data_manager.get_segments(we_iter.n_iter)
        elif opts.traj_type == 'complete':
            segments = data_manager.get_segments(1, n_iter_upper = we_iter.n_iter, 
                                                 status_criteria = Segment.SEG_STATUS_COMPLETE)
        elif opts.traj_type == 'merged':
            segments = data_manager.get_segments(1, n_iter_upper = we_iter.n_iter, 
                                                 endpoint_criteria = Segment.SEG_ENDPOINT_TYPE_MERGED)    
        elif opts.traj_type == 'recycled':
            segments = data_manager.get_segments(1, n_iter_upper = we_iter.n_iter, 
                                                 endpoint_criteria = Segment.SEG_ENDPOINT_TYPE_RECYCLED)             
        elif opts.traj_type == 'all':
            segments = data_manager.get_segments(1, n_iter_upper = we_iter.n_iter)

        for segment in segments:
            sys.stdout.write('%-12d     %-12d    %21.16g\n'
                                     % (segment.seg_id, segment.n_iter,
                                        segment.weight))
                
    def cmd_transanl(self, args):
        # Find transitions
        parser = self.make_parser()
        self.add_iter_param(parser)
        parser.add_option('-a', '--analysis-config', dest='anl_config',
                          help='use ANALYSIS_CONFIG as configuration file '
                              +'(default: analysis.cfg)')
        parser.add_option('-o', '--output', dest='output_pattern',
                          help='write results to OUTPUT_PATTERN, which must '
                              +'contain two "%s" flags which will be replaced '
                              +'by the region names as given in the analysis '
                              +'configuration file')
        parser.set_defaults(anl_config = 'analysis.cfg',
                            output_pattern = 'ed_%s_%s.txt')
        (opts,args) = parser.parse_args(args)
        
        from wemd.util.config_dict import ConfigDict, ConfigError
        from wemd.analysis.transitions import TransitionEventAccumulator, OneDimRegionSet
        from wemd.analysis.trajtree import TrajTree

        transcfg = ConfigDict({'data.squeeze': False,
                               'data.timestep.units': 'ps'})
        transcfg.read_config_file(opts.anl_config)
        transcfg.require('regions.edges')
        transcfg.require('regions.names')
        
        region_names = transcfg.get_list('regions.names')
        region_edges = [float(v) for v in transcfg.get_list('regions.edges')]
        if len(region_edges) != len(region_names) + 1:
            sys.stderr.write('region names and boundaries do not match\n')
            sys.exit(EX_ERROR)
            
        try:
            translog = transcfg.get_file_object('output.transition_log', mode='w')
        except KeyError:
            translog = None 
        
        squeeze_data = transcfg.get_bool('data.squeeze')
        region_boundaries = []
        for (irr,rname) in enumerate(region_names):
            region_boundaries.append((region_edges[irr], region_edges[irr+1]))
        regions = OneDimRegionSet(region_names, region_boundaries)
        
        self.load_sim_manager()
        self.sim_manager.load_data_manager()
        
        final_we_iter = self.get_sim_iter(opts.we_iter, complete = True)
        max_iter = final_we_iter.n_iter
        data_manager = self.sim_manager.data_manager

        tree = TrajTree(data_manager, False)
        trans_finder = TransitionEventAccumulator(regions, data_overlaps = (not squeeze_data))
        
        def analyze_segment(segment, children, history):
            seg_weights = numpy.zeros((len(segment.pcoord),),numpy.float64)
            seg_weights += segment.weight
            trans_finder.timestep = segment.data.get('dt', 1.0)
            trans_finder.identify_transitions(segment.pcoord, seg_weights)
            
        import numpy
        
        tree.trace_trajectories(max_iter, analyze_segment, 
                                trans_finder.get_state, 
                                trans_finder.set_state)
        sys.stdout.write('event count (row->column, states %s\n' % ', '.join(regions.names))
        sys.stdout.write('%s\n' % trans_finder.event_counts)
        for ((ir1, ir2), ed_list) in trans_finder.eds.iteritems():
            region1_name = regions.names[ir1]
            region2_name = regions.names[ir2]
            if len(ed_list) == 0:
                sys.stdout.write('No %s->%s transitions observed\n'
                                         % (region1_name, region2_name))
            else:
                ed_array = numpy.array(ed_list, numpy.float64)
                ed_array[:,0] *= trans_finder.timestep
                sys.stdout.write('\nStatistics for %s->%s:\n'
                                         % (region1_name, region2_name))
                (ed_mean, ed_norm) = numpy.average(ed_array[:,0],
                                                   weights = ed_array[:,1],
                                                   returned = True)    
                sys.stdout.write('Number of events:    %d\n' % ed_array.shape[0])
                sys.stdout.write('ED average:          %g\n' % ed_array[:,0].mean())
                sys.stdout.write('ED weighted average: %g\n' % ed_mean)
                sys.stdout.write('ED min:              %g\n' % ed_array[:,0].min())
                sys.stdout.write('ED median:           %g\n' % numpy.median(ed_array[:,0]))
                sys.stdout.write('ED max:              %g\n' % ed_array[:,0].max())
                
                ed_file = open(opts.output_pattern % 
                               (region1_name, region2_name), 'wt')
                for irow in xrange(0, ed_array.shape[0]):
                    ed_file.write('%20.16g    %20.16g\n'
                                  % tuple(ed_array[irow,0:2]))
                ed_file.close()
                    
            
    def cmd_tracetraj(self, args):
        parser = self.make_parser('TRAJ_ID', description = 'trace trajectory '
                                  +'path')
        (opts, args) = parser.parse_args(args)
        if len(args) != 1:
            parser.print_help(sys.stderr)
            sys.exit(EX_USAGE_ERROR)
        
        self.load_sim_manager()
        self.sim_manager.load_data_manager()
        data_manager = self.sim_manager.data_manager
        seg = data_manager.get_segment(long(args[0]), load_p_parent = True)
        
        segs = [seg]
        while seg.p_parent and seg.n_iter:
            seg = data_manager.get_segment(seg.p_parent.seg_id, load_p_parent = True)
            segs.append(seg)
        traj = list(reversed(segs))
        
        for seg in traj:
            
            if seg.p_parent is not None:
                p_parent_id = seg.p_parent.seg_id
            else:
                p_parent_id = 0
            
            if seg.pcoord is not None:
                pcoord = seg.pcoord[-1]
            else:
                pcoord = -1 #ie an incomplete seg was given
                
            sys.stdout.write('%5d  %10d  %10d  %21.16g  %21.16g\n'
                                     % (seg.n_iter,
                                        seg.seg_id,
                                        p_parent_id,
                                        seg.weight,
                                        pcoord))
    def cmd_fluxanl(self, args):
        parser = self.make_parser()
        parser.add_option('-T', '--tau', dest='tau', type='float',
                          default=1.0,
                          help='length of each WE iteration in simulation '
                              +'time is TAU (default: 1.0)')
        (opts,args) = parser.parse_args(args)
        self.load_sim_manager()
        self.sim_manager.load_data_manager()
        data_manager = self.sim_manager.data_manager
        
        latest_we_iter = self.get_sim_iter(None, complete = True)

        import numpy
        fluxen = numpy.zeros((latest_we_iter.n_iter,), numpy.float64)
        for i_iter in xrange(1, latest_we_iter.n_iter+1):
            we_iter = self.get_sim_iter(i_iter)
            fluxen[i_iter-1] = we_iter.data['recycled_population']
            
        fluxen /= opts.tau
        
        for irow in xrange(0, fluxen.shape[0]):
            sys.stdout.write('%-8d    %16.12g    %21.16g\n'
                                     % (irow+1, opts.tau*(irow+1), 
                                        fluxen[irow]))
        
            
            
            
        
        
        
if __name__ == '__main__':
    WEMDAnlTool().run()

