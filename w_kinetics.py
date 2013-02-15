'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from west.data_manager import seg_id_dtype
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection, BinMappingComponent
import numpy, h5py
from westtools import h5io

log = logging.getLogger('westtools.w_kinetics')

class KineticsCrawler:
    def __init__(self):
        pass

def assign_trajectories(data_reader, mapper, output_group, output_name='bin_assignments'):
    '''Assign all trajectories in the WEST HDF5 file handled by ``data_reader``
    to bins using the given ``mapper``. Results are stored in ``output_group``
    as an integer dataset ``bin_assignments[n_iter][seg_id][timepoint]``.'''
    
    iter_count = data_reader.current_iteration
    max_seg_count = 0
    max_timepoints = 0
    
    assignments_dtype = numpy.min_scalar_type(mapper.nbins)
    
    # Determine dimensions of output array
    for n_iter in xrange(1, iter_count):
        iter_group = data_reader.get_iter_group(n_iter)
        pcoord_shape = iter_group['pcoord'].shape
        max_seg_count = max(max_seg_count, pcoord_shape[0])
        max_timepoints = max(max_timepoints, pcoord_shape[1])
        del iter_group
            
    assignments_shape = (iter_count,max_seg_count,max_timepoints)
    
    try:
        del output_group[output_name]
    except KeyError:
        pass
    
    # Create output array
    assignments_ds = output_group.create_dataset(output_name,
                                                 shape=assignments_shape,
                                                 chunks=h5io.calc_chunksize(assignments_shape, assignments_dtype),
                                                 dtype=assignments_dtype,
                                                 shuffle=True,
                                                 compression=9)
    
    # Assign trajectories
    for (iiter, n_iter) in enumerate(xrange(1,iter_count)):
        print('iteration {:d}'.format(n_iter))
        iter_group = data_reader.get_iter_group(n_iter)
        pcoords = iter_group['pcoord'][...]
        assignments = numpy.empty(pcoords.shape[:-1], assignments_dtype)
        for iseg in xrange(pcoords.shape[0]):
            assignments[iseg,:] = mapper.assign(pcoords[iseg])
        assignments_ds[iiter,:pcoords.shape[0],...] = assignments
        del pcoords, iter_group, assignments
    

class WKinetics(WESTTool):
    prog='w_kinetics'
    description = '''\
Calculate kinetics data for arbitrary states from weighted ensemble data.
Bins here represent kinetic macrostates (e.g. initial and final states),
not microstates (e.g. fine-grained bins used during a WE simulation). 
'''
    
    def __init__(self):
        super(WKinetics,self).__init__()
        self.data_reader = WESTDataReader() 
        self.binning = BinMappingComponent()
        self.output_file = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.binning.add_args(parser, suppress=['--bins-from-system'])
        
        parser.add_argument('-o', '--output', dest='output', default='kinetics.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        
    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open(mode='r+')
        self.binning.process_args(args)
        self.output_file = h5py.File(args.output, 'w')
        h5io.stamp_creator_data(self.output_file)
        
    def go(self):
        assign_trajectories(self.data_reader, self.binning.mapper, self.output_file)
        

if __name__ == '__main__':
    WKinetics().main()
