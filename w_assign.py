from __future__ import print_function, division; __metaclass__ = type
import logging
import math
from west.data_manager import seg_id_dtype
from west.binning.assign import index_dtype
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection, BinMappingComponent
import numpy, h5py
from westtools import h5io
from westtools.h5io import WESTPAH5File
from westpa.extloader import get_object

log = logging.getLogger('westtools.w_assign')

def parse_pcoord_value(pc_str):
    namespace = {'math': math,
                 'numpy': numpy,
                 'inf': float('inf')}
    
    arr = numpy.array(eval(pc_str,namespace))
    if arr.ndim == 0:
        arr.shape = (1,1)
    elif arr.ndim == 1:
        arr.shape = (1,) + arr.shape 
    else:
        raise ValueError('too many dimensions')
    return arr


def default_construct_pcoord(n_iter, iter_group):
    return iter_group['pcoord'][...]
    

class WAssign(WESTTool):
    prog='w_assign'
    description = '''\
Assign walkers to bins, producing a file (by default named "assign.h5")
which can be used in subsequent analysis.

Progress coordinate data is taken by default from the "pcoord" dataset
for each iteration in the main HDF5 file (usually west.h5). However,
an arbitrary function can be provided with -p/--construct-pcoord,
which can be used to consolidate data from several sources or otherwise
preprocess it before binning occurs. 

Optionally, a list of coordinate tuples may be provided. These tuples 
will be mapped to bins, and those bins recorded as kinetic macrostates
to be used in subsequent kinetics analysis. 

'''
    
    def __init__(self):
        super(WAssign,self).__init__()
        self.data_reader = WESTDataReader() 
        self.binning = BinMappingComponent()
        self.iter_range = IterRangeSelection()
        self.output_file = None
        self.construct_pcoord = default_construct_pcoord
        self.state_points = []
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        self.binning.add_args(parser, suppress=['--bins-from-file'])
        self.iter_range.add_args(parser)
        
        agroup = parser.add_argument_group('other options')
        
        agroup.add_argument('-p', '--construct-pcoord', metavar='MODULE.FUNCTION',
                             help='''Use the given function to construct progress coordinate data
                             for each iteration. This function will be called once per iteration as
                             ``get_pcoord(n_iter, iter_group)``, and must return an array indexable as
                             [seg_id][timepoint][dimension]. By default, the "pcoord" dataset is
                             loaded into memory and returned.''')
        
        agroup.add_argument('-o', '--output', dest='output', default='assign.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        
        agroup.add_argument('state_points', metavar='PC_TUPLE',
                            nargs='*',
                            help='Coordinates of kinetic macrostates (optional).')

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open(mode='r')
        self.iter_range.process_args(args)
        self.binning.process_args(args)
        if args.construct_pcoord:
            self.construct_pcoord = get_object(args.construct_pcoord,path=['.'])
        if args.state_points:
            self.state_points.extend([parse_pcoord_value(pcstr) for pcstr in args.state_points])
        self.output_file = WESTPAH5File(args.output, 'w', creating_program=True)
        
    def go(self):
        assign = self.binning.mapper.assign
        
        iter_start = self.iter_range.iter_start 
        iter_stop =  self.iter_range.iter_stop
        
        h5io.stamp_iter_range(self.output_file, iter_start, iter_stop)
        
        self.output_file.attrs['nbins'] = self.binning.mapper.nbins
        
        if self.state_points:
            state_points_array = numpy.vstack(self.state_points)
            state_assignments_array = assign(state_points_array)
            
            self.output_file['state_points'] = state_points_array
            self.output_file['state_assignments'] = state_assignments_array 
        
        for n_iter in xrange(iter_start,iter_stop):
            iter_group = self.data_reader.get_iter_group(n_iter)
            pcoords = self.construct_pcoord(n_iter,iter_group)
            
            n_segs, n_pts, n_dim = pcoords.shape
            
            output_iter_group = self.output_file.create_iter_group(n_iter)
            
            assignments = numpy.empty((n_segs,n_pts), index_dtype)
            mask = numpy.ones((n_pts,), numpy.bool_)
            
            for seg_id in xrange(n_segs):
                assign(pcoords[seg_id], mask, assignments[seg_id])
                
            output_iter_group.create_dataset('assignments',
                                             data=assignments, shuffle=True, compression=9,
                                             chunks = h5io.calc_chunksize(assignments.shape, assignments.dtype))
            
            del pcoords, assignments, mask, iter_group

if __name__ == '__main__':
    WAssign().main()

