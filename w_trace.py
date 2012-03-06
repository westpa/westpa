from __future__ import print_function, division; __metaclass__ = type
import sys
from wt2.tool_classes import WEMDTool, HDF5Storage, WEMDDataReader
import wemd
import numpy

from wemd import Segment
from wemd.states import InitialState
from wemd.data_manager import weight_dtype, n_iter_dtype, seg_id_dtype, utime_dtype

class Trace:
    '''A class representing a trace of a certain trajectory segment back to its origin.'''
    
    def __init__(self, summary, endpoint_type, basis_state, initial_state, data_manager = None):
        self.summary = summary
        self.endpoint_type = endpoint_type
        self.basis_state = basis_state
        self.initial_state = initial_state
        self.data_manager = data_manager or wemd.rc.get_data_manager()
            
    def __len__(self):
        try:
            return len(self.summary)
        except TypeError:
            return 0
        
    def __getitem__(self, sl):
        return self.summary[sl]
    
    def __iter__(self):
        return iter(self.summary)
            
    @classmethod
    def from_data_manager(cls, n_iter, seg_id, data_manager = None):
        '''Construct and return a trajectory trace whose last segment is identified
        by ``seg_id`` in the iteration number ``n_iter``.'''
        
        data_manager = data_manager or wemd.rc.get_data_manager()
        
        # These values are used later on
        endpoint_type = None
        pcoord_dtype = None
        pcoord_pt_shape = None
        
        seginfo = []
        parent_id = seg_id
        
        while n_iter > 0 and parent_id >= 0:
            seg_id = parent_id
            iter_group = data_manager.get_iter_group(n_iter)
            pcoord_ds = iter_group['pcoord']
            seg_index = iter_group['seg_index']
            n_segs = pcoord_ds.shape[0]
            pcoord_len = pcoord_ds.shape[1]
            
            assert seg_id < n_segs
            
            indexrow = seg_index[seg_id]
            final_pcoord = pcoord_ds[seg_id, pcoord_len-1]
            weight = indexrow['weight']
            cputime = indexrow['cputime']
            walltime = indexrow['walltime']
            
            try:
                parent_id = long(indexrow['parent_id'])
            except IndexError:
                # old HDF5 version
                parent_id = long(iter_group['parents'][indexrow['parents_offset']])
                
            if endpoint_type is None:
                endpoint_type = indexrow['endpoint_type']
                pcoord_pt_shape = pcoord_ds.shape[2:]
                pcoord_dtype = pcoord_ds.dtype
                
            seginfo.append((n_iter, seg_id, weight, walltime, cputime, final_pcoord))
            
            del iter_group, pcoord_ds, seg_index
            n_iter -= 1
            
        # loop terminates with parent_id set to the identifier of the initial state, 
        # seg_id set to the identifier of the first segment in the trajectory, and
        # n_iter set to one less than the iteration of the first segment
        first_iter = n_iter + 1
        first_seg_id = seg_id
        first_parent_id = parent_id
        
        seginfo.reverse()
        
        summary_dtype = numpy.dtype([('n_iter', n_iter_dtype),
                                     ('seg_id', seg_id_dtype),
                                     ('weight', weight_dtype),
                                     ('walltime', utime_dtype),
                                     ('cputime', utime_dtype),
                                     ('final_pcoord', pcoord_dtype, pcoord_pt_shape),
                                     ])
        
        summary = numpy.array(seginfo, dtype=summary_dtype)
        
        try:
            initial_state = data_manager.get_segment_initial_states([first_seg_id], first_iter)[0]
        except KeyError:
            # old HDF5 version
            assert parent_id < 0
            istate_pcoord = data_manager.get_iter_group(first_iter)['pcoord'][first_seg_id,0]
            istate_id = -(first_parent_id+1)
            basis_state = None
            initial_state = InitialState(istate_id, None, iter_created=0, pcoord=istate_pcoord)
            
        else:
            basis_state = data_manager.get_basis_states(first_iter)[initial_state.basis_state_id]
            
        return cls(summary, endpoint_type, basis_state, initial_state, data_manager)

    def trace_timepoint_dataset(self, dsname):
        '''Return a trace along this trajectory over a dataset which is layed out as [seg_id][timepoint][...].
        Overlapping values at segment boundaries are accounted for.  Returns (data_trace, weight), where 
        data_trace is a time series of the dataset along this trajectory, and weight is the corresponding
        trajectory weight at each time point.'''
        
        first_n_iter, first_seg_id = self.summary[0]['n_iter'], self.summary[0]['seg_id']
        first_iter_group = self.data_manager.get_iter_group(first_n_iter)
        first_iter_ds = first_iter_group[dsname]
        n_segs = len(self)
        n_points_per_seg = first_iter_ds.shape[1]
        
        length = n_points_per_seg + (n_segs-1)*(n_points_per_seg-1)
        print(length) 
        tracedata = numpy.empty((length,) + first_iter_ds.shape[2:], dtype=first_iter_ds.dtype)
        traceweight = numpy.empty((length,), weight_dtype)
        
        # Store first segment data
        tracedata[0:n_points_per_seg] = first_iter_ds[first_seg_id]
        traceweight[0:n_points_per_seg] = self.summary[0]['weight']
        
        # Store remainder of data
        for iseg, summary_item in enumerate(self.summary[1:]):
            n_iter = summary_item['n_iter']
            seg_id = summary_item['seg_id']
            iter_group = self.data_manager.get_iter_group(n_iter)
            offset = n_points_per_seg + iseg*(n_points_per_seg-1)
            length = n_points_per_seg - 1
            seg_data = iter_group[dsname][seg_id][1:]
            weight = summary_item['weight']
            tracedata[offset:offset+length] = seg_data 
            traceweight[offset:offset+length] = weight
            del seg_data, iter_group
            
        return tracedata, traceweight
    
    def trace_perseg_dataset(self, dsname):
        '''Return a trace along this trajectory over a dataset which is layed out as [seg_id][...].
        Returns (data_trace, weight), where  data_trace is a time series of the dataset along this
        trajectory, and weight is the corresponding trajectory weight at each time point.'''

        first_n_iter, first_seg_id = self.summary[0]['n_iter'], self.summary[0]['seg_id']
        first_iter_group = self.data_manager.get_iter_group(first_n_iter)
        first_iter_ds = first_iter_group[dsname]
        n_segs = len(self)
        tracedata = numpy.empty((n_segs,) + first_iter_ds.shape[1:], dtype=first_iter_ds.dtype)
        traceweight = numpy.empty((n_segs,), weight_dtype)
        tracedata[0] = first_iter_ds[first_seg_id]
        traceweight[0] = self.summary[0]['weight']
        for isegm1, summary_item in enumerate(self.summary[1:]):
            iseg = isegm1 + 1
            n_iter = summary_item['n_iter']
            seg_id = summary_item['seg_id']            
            iter_group = self.data_manager.get_iter_group(n_iter)
            seg_data = iter_group[dsname][seg_id]
            tracedata[iseg] = seg_data
            traceweight[iseg] = summary_item['weight']
            del seg_data
            
        return tracedata, traceweight
        
class WTraceTool(WEMDTool):
    prog='w_trace'
    description = '''\
Trace individual WEMD trajectories and emit (or calculate) quantities along the
trajectory.

Trajectories are specified as N_ITER:SEG_ID pairs. Each segment is traced back
to its initial point, and then various quantities (notably n_iter and seg_id)
are printed in order from initial point up until the given segment in the given
iteration.

Output is stored in several files, all named according to the pattern given by
the -o/--output-pattern parameter. The default output pattern is "traj_%d_%d",
where the printf-style format codes are replaced by the iteration number and
segment ID of the terminal segment of the trajectory being traced.
'''


    pcoord_formats = {'u8': '%20d',
                      'i8': '%20d',
                      'u4': '%10d',
                      'i4': '%11d',
                      'u2': '%5d',
                      'i2': '%6d',
                      'f4': '%14.7g',
                      'f8': '%023.15g'}
    
    def __init__(self):
        super(WTraceTool,self).__init__()
        
        self.data_reader = WEMDDataReader()
        self.h5storage = HDF5Storage()
        self.output_pattern = None
        self.endpoints = None

        
    # Interface for command-line tools
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.h5storage.add_args(parser)
        parser.add_argument('endpoints',  metavar='N_ITER:SEG_ID', nargs='+',
                            help='''Trace trajectory ending (or at least alive at) N_ITER:SEG_ID.''')
        
        #tgroup = parser.add_argument_group('trace options')
        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output-pattern', default='traj_%d_%d',
                            help='''Write per-trajectory data to output files whose names begin with OUTPUT_PATTERN, which
                                 must contain two printf-style format flags which will be replaced with the iteration number
                                 and segment ID of the terminal segment of the trajectory being traced.
                                 (Default: %(default)s.)''')
        
    
    def process_args(self, args):
        self.data_reader.process_args(args)
        self.h5storage.process_args(args)
        self.endpoints = [map(long,endpoint.split(':')) for endpoint in args.endpoints]
        self.output_pattern = args.output_pattern
        
        self.h5storage.open_analysis_h5file()
    
    def go(self):
        self.data_reader.open('r')
        
        for n_iter, seg_id in self.endpoints:
            trace = Trace.from_data_manager(n_iter,seg_id, self.data_reader.data_manager)
            #tsum = tracer.summarize_trace(trace)
            #segs = tracer.segments_from_trace(trace,load_pcoords=True)
            #print(trace.initial_state)
            #print(trace.summary)
            #print(trace.trace_timepoint_dataset('pcoord'))
            #with open(self.output_pattern % (n_iter,seg_id) + '_trace.txt', 'wt') as trace_output:
            #    self.emit_trace_text(segs, trace_output)
            self.emit_trace_text(trace, sys.stdout)
    
    def emit_trace_h5(self, traced_segs, output_group, trace_ds_name='trace'):
        '''Dump summary information about each segment in the given trace to the given group in the 
        current analysis HDF5 file.  A table named ``trace_ds_name`` is created.'''
        pass
        
            
    def emit_trace_text(self, trace, output_file):
        '''Dump summary information about each segment in the given trace to the given output_file,
        which must be opened for writing in text mode.  Output columns are separated by at least
        one space.'''
        
        if not trace:
            return
                
        pcoord_ndim = trace[0]['final_pcoord'].shape[0]
        lastseg = trace[-1]
        len_n_iter = max(6, len(str(lastseg['n_iter'])))
        len_seg_id = max(6, max(len(str(seg_id)) for seg_id in trace['seg_id']))
        seg_pattern = '    '.join(['{n_iter:{len_n_iter}d}',
                                   '{seg_id:{len_seg_id}d}',
                                   '{weight:22.17e}',
                                   '{walltime:10.6g}',
                                   '{cputime:10.6g}',
                                   '{pcoord_str:s}'
                                   ]) + '\n'
                                  
        
        output_file.write('''\
# Trace of trajectory ending in n_iter:seg_id {n_iter:d}:{seg_id:d} (endpoint type {endpoint_type_text:s})   
# column  0: iteration (0 => initial state)
# column  1: seg_id
# column  2: weight
# column  3: wallclock time (s)
# column  4: CPU time (s)
'''.format(n_iter = long(lastseg['n_iter']), 
           seg_id = long(lastseg['seg_id']), 
           endpoint_type_text = Segment.endpoint_type_names[trace.endpoint_type]))
        
        
        if pcoord_ndim == 1:
            output_file.write('''\
# column  5: final progress coordinate value            
''')
        else:
            fpcbegin = 5
            fpcend = fpcbegin + pcoord_ndim - 1
            output_file.write('''\
# columns {fpcbegin:d} -- {fpcend:d}: final progress coordinate value
'''.format(fpcbegin=fpcbegin,fpcend=fpcend))
        
        
        pcoord_formats = self.pcoord_formats
        for segment in trace:            
            pcoord_str = '    '.join(pcoord_formats.get(pcfield.dtype.str[1:], '%s') % pcfield 
                                     for pcfield in segment['final_pcoord'])
            output_file.write(seg_pattern.format(n_iter = long(segment['n_iter']), 
                                                 seg_id = long(segment['seg_id']),
                                                 weight = float(segment['weight']),
                                                 walltime = float(segment['walltime']),
                                                 cputime = float(segment['cputime']), 
                                                 pcoord_str=pcoord_str,
                                                 len_n_iter=len_n_iter,
                                                 len_seg_id=len_seg_id))
        
    
        

if __name__ == '__main__':
    WTraceTool().main()
    
