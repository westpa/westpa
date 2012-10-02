from __future__ import division, print_function
import sys, argparse
import numpy
import west, oldtools

from west import Segment
from oldtools.aframe import WESTAnalysisTool, WESTDataReaderMixin, CommonOutputMixin

import logging
log = logging.getLogger('w_trace')

class WTrace(CommonOutputMixin,WESTDataReaderMixin,WESTAnalysisTool):
    def __init__(self):
        super(WTrace,self).__init__()
        self.output_pattern = None
        self.terminal_segments = []
    
    def trace_trajs(self):
        
        report_format = '    '.join(['{segment.n_iter:8d}',
                                     '{segment.seg_id:8d}',
                                     '{segment.weight:20.14g}',
                                     ])
        
        # DAMMIT NUMPY GET YOUR ACT TOGETHER WITH STRING FORMATTING
#        pcoord_formats = {'u8': '{:20d}',
#                          'i8': '{:20d}',
#                          'u4': '{:10d}',
#                          'i4': '{:11d}',
#                          'u2': '{:5d}',
#                          'i2': '{:6d}',
#                          'f4': '{:14.7g}',
#                          'f8': '{:23.15g}'}

        pcoord_formats = {'u8': '%20d',
                          'i8': '%20d',
                          'u4': '%10d',
                          'i4': '%11d',
                          'u2': '%5d',
                          'i2': '%6d',
                          'f4': '%14.7g',
                          'f8': '%23.15g'}

                         
        for (n_iter,traj_seg_id) in self.terminal_segments:
            output_file = file(self.output_pattern.format(n_iter=n_iter,seg_id=traj_seg_id),'wt')
            if not self.output_suppress_headers:
                output_file.write('''\
# Trace of trajectory terminated by segment {n_iter:d}:{traj_seg_id:d}
# column 0:    n_iter (0 represents initial position)
# column 1:    seg_id
# column 2:    weight
# column>2:    final progress coordinate value in given segment
'''.format(n_iter=n_iter,traj_seg_id=traj_seg_id))


            segments = []
            seg_id = traj_seg_id
            while seg_id >= 0:                
                iter_group = self.get_iter_group(n_iter)
                seg_info = iter_group['seg_index'][seg_id] 
                final_pcoord = iter_group['pcoord'][seg_id,-1,:]
                p_parent_id = seg_info['parent_id']                 
                segments.append(Segment(n_iter=n_iter, seg_id=seg_id, 
                                        weight=seg_info['weight'],
                                        parent_id=p_parent_id,
                                        pcoord=final_pcoord))                
                seg_id = p_parent_id
                n_iter -= 1
            
            # Add a fictitious segment for the initial position
            segments.append(Segment(n_iter=0, seg_id=seg_id, weight=segments[-1].weight,
                                    pcoord=iter_group['pcoord'][segments[-1].seg_id,0,:]))    
            
            # for segment in segments: print
            for segment in reversed(segments):
                output_file.write(report_format.format(segment=segment))
                fields = ['']
                for field in segment.pcoord:
                    fields.append(pcoord_formats.get(field.dtype.str[1:], '%s') % field)
                output_file.write('    '.join(fields))
                output_file.write('\n')

wtrace = WTrace()

parser = argparse.ArgumentParser('w_trace', description='''\
Trace trajectories. One or more trajectories must be specified, each as n_iter:seg_id''')
west.rc.add_args(parser)
wtrace.add_args(parser)
parser.add_argument('-o', '--output', dest='output_pattern', default='traj_{n_iter:06d}_{seg_id:06d}.txt',
                    help='''Store output in OUTPUT_PATTERN, which must be a Python3-style format
                         containing the named fields 'n_iter' and 'seg_id' (default: '%(default)s').''')
parser.add_argument('segments', nargs='+', metavar='SEGMENT',
                    help='Segment(s) to trace, each specified as "n_iter:seg_id"')
args = parser.parse_args()
west.rc.process_args(args, config_required=False)
wtrace.process_args(args)

wtrace.output_pattern = args.output_pattern

for segspec in args.segments:
    try:
        n_iter,seg_id = map(long,segspec.split(':'))
    except Exception as e:
        sys.stderr.write('Invalid terminal segment specification {!r}\n'.format(segspec))
        sys.exit(3)
    else:
        wtrace.terminal_segments.append((n_iter,seg_id))
wtrace.trace_trajs()

#ttree = oldtools.trajectories.trajtree.TrajTree(data_manager)
#for segspec in args.segment:
#    n_iter, seg_id = map(int, segspec.split(':'))
#    trajectory = ttree.trace_to_root(n_iter, seg_id)
#    osegspec = '{n_iter:d}:{seg_id:d}'.format(n_iter = n_iter, seg_id = seg_id)
#    for segment in trajectory:
#        output_file.write('{osegspec:20s}    {segment.n_iter:10d}    {segment.seg_id:10d}\n'
#                          .format(osegspec = osegspec, segment = segment))     
