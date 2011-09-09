from __future__ import print_function, division; __metaclass__=type

import logging
log = logging.getLogger(__name__)

import wemd, wemdtools
from wemdtools.aframe import AnalysisMixin, ArgumentError
import numpy
from wemd import Segment

class DataReaderMixin(AnalysisMixin):
    '''A mixin for analysis requiring access to the HDF5 files generated during a WEMD run.'''

    def __init__(self):
        super(DataReaderMixin,self).__init__()
        
        self.data_manager = None
        self.run_h5name = None
        
        # Whether pcoord caching is active
        self.__cache_pcoords = False        

        # Cached items
        self.__c_summary = None
        self.__c_iter_groups = dict()
        self.__c_seg_id_ranges = dict()
        self.__c_seg_indices = dict()
        self.__c_parent_arrays = dict()
        self.__c_p_parent_arrays = dict()
        self.__c_pcoord_arrays = dict()
        self.__c_pcoord_datasets = dict()
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upcall = super(DataReaderMixin,self).add_args
            except AttributeError:
                pass
            else:
                upcall(parser)
        
        group = parser.add_argument_group('WEMD HDF5 options')
        #group.add_argument('--no-cache', dest='no_cache_rundata', action='store_true',
        #                    help='''Disable all caching of WEMD run data.''')
        group.add_argument('run_h5name', nargs='?', metavar='WEMD_H5FILE',
                           help='''Take data from WEMD_H5FILE (default: read from the HDF5 file specified in wemd.cfg).''')

    def process_args(self, args, upcall = True):            
        if args.run_h5name:
            self.run_h5name = args.run_h5name
        else:
            wemd.rc.config.require('data.h5file')
            self.run_h5name = wemd.rc.config.get_path('data.h5file') 
        
        wemd.rc.pstatus("Using run data from '{}'".format(self.run_h5name))
        
        self.data_manager = wemd.rc.get_data_manager()
        self.data_manager.backing_file = self.run_h5name        
        self.data_manager.open_backing(mode='r')
        
        if upcall:
            try:
                upfunc = super(DataReaderMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)

    def clear_run_cache(self):
        del self.__c_summary
        del self.__c_iter_groups, self.__c_seg_id_ranges, self.__c_seg_indices, self.__c_parent_arrays, self.__c_p_parent_arrays
        del self.__c_pcoord_arrays, self.__c_pcoord_datasets
        
        self.__c_summary = None
        self.__c_iter_groups = dict()
        self.__c_seg_id_ranges = dict()
        self.__c_seg_indices = dict()
        self.__c_parent_arrays = dict()
        self.__c_p_parent_arrays = dict()
        self.__c_pcoord_arrays = dict()
        self.__c_pcoord_datasets = dict()

    @property
    def cache_pcoords(self):
        '''Whether or not to cache progress coordinate data. While caching this data
        can significantly speed up some analysis operations, this requires
        copious RAM.
        
        Setting this to False when it was formerly True will release any cached data. 
        '''
        return self.__cache_pcoords

    @cache_pcoords.setter
    def cache_pcoords(self, cache):
        self.__cache_pcoords = cache
        
        if not cache:
            del self.__c_pcoord_arrays
            self.__c_pcoord_arrays = dict()
            
    def get_summary_table(self):
        if self.__c_summary is None:
            self.__c_summary = self.data_manager.h5file['/summary'][...]
        return self.__c_summary

    def get_iter_group(self, n_iter):
        '''Return the HDF5 group corresponding to ``n_iter``'''
        try:
            return self.__c_iter_groups[n_iter]
        except KeyError:
            iter_group = self.data_manager.get_iter_group(n_iter)
            return iter_group
    
    def get_segments(self, n_iter, include_pcoords = True):
        '''Return all segments present in iteration n_iter'''
        return self.get_segments_by_id(n_iter, self.get_seg_ids(n_iter, None), include_pcoords)
    
    def get_segments_by_id(self, n_iter, seg_ids, include_pcoords = True):
        '''Get segments from the data manager, employing caching where possible'''
        
        if len(seg_ids) == 0: return []
        
        seg_index  = self.get_seg_index(n_iter)
        all_parent_ids = self.get_parent_array(n_iter)
        
        segments = []

        if include_pcoords:
            pcoords = self.get_pcoords(n_iter, seg_ids)
        
        for (isegid, seg_id) in enumerate(seg_ids):
            row = seg_index[seg_id]
            parents_offset = row['parents_offset']
            n_parents = row['n_parents']
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              n_parents = n_parents,
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'],
                              )
            if include_pcoords:
                segment.pcoord = pcoords[isegid]

            parent_ids = all_parent_ids[parents_offset:parents_offset+n_parents]
            segment.parent_ids = {long(parent_id) for parent_id in parent_ids}
            segment.p_parent_id = long(parent_ids[0])
            segments.append(segment)

        return segments        
    
    def get_children(self, segment, include_pcoords=True):
        p_parents = self.get_p_parent_array(segment.n_iter+1)
        seg_ids = self.get_seg_ids(segment.n_iter+1, p_parents == segment.seg_id)
        return self.get_segments_by_id(segment.n_iter+1, seg_ids, include_pcoords)

    def get_seg_index(self, n_iter):
        try:
            return self.__c_seg_indices[n_iter]
        except KeyError:
            seg_index = self.__c_seg_indices[n_iter] = self.get_iter_group(n_iter)['seg_index'][...]
            return seg_index
        
    def get_parent_array(self, n_iter):
        try:
            return self.__c_parent_arrays[n_iter]
        except KeyError:
            parent_array = self.__c_parent_arrays[n_iter] = self.get_iter_group(n_iter)['parents'][...]
            return parent_array
                
    def get_p_parent_array(self, n_iter):
        try:
            return self.__c_p_parent_arrays[n_iter]
        except KeyError:
            parent_offsets = self.get_seg_index(n_iter)['parents_offset'][:]
            p_parent_array = self.get_parent_array(n_iter)[parent_offsets]
            self.__c_p_parent_arrays[n_iter] = p_parent_array
            return p_parent_array
                
    def get_pcoord_array(self, n_iter):
        try:
            return self.__c_pcoord_arrays[n_iter]
        except KeyError:
            pcoords = self.__c_pcoord_arrays[n_iter] = self.get_iter_group(n_iter)['pcoord'][...]
            return pcoords
        
    def get_pcoord_dataset(self, n_iter):
        try:
            return self.__c_pcoord_datasets[n_iter]
        except KeyError:
            pcoord_ds = self.__c_pcoord_datasets[n_iter] = self.get_iter_group(n_iter)['pcoord']
            return pcoord_ds
        
    def get_pcoords(self, n_iter, seg_ids):
        if self.__cache_pcoords:
            pcarray = self.get_pcoord_array(n_iter)
            return [pcarray[seg_id,...] for seg_id in seg_ids]
        else:
            return self.get_pcoord_dataset(n_iter)[list(seg_ids),...]
        
    def get_seg_ids(self, n_iter, bool_array = None):
        try:
            all_ids = self.__c_seg_id_ranges[n_iter]
        except KeyError:
            all_ids = self.__c_seg_id_ranges[n_iter] = numpy.arange(0,len(self.get_seg_index(n_iter)), dtype=numpy.uint32)
            
            
        if bool_array is None:             
            return all_ids
        else:        
            seg_ids = all_ids[bool_array]        
            try:
                if len(seg_ids) == 0: return []
            except TypeError:
                # Not iterable
                return [seg_ids]
            else:
                return seg_ids
    
    def get_created_seg_ids(self, n_iter):
        '''Return a list of seg_ids corresponding to segments which were created for the given iteration (are not
        continuations).'''
        
        # Created segments have p_parent_id < 0
        p_parent_ids = self.get_p_parent_array(n_iter)        
        return self.get_seg_ids(n_iter, p_parent_ids < 0)


    def max_iter_segs_in_range(self, first_iter, last_iter):
        '''Return the maximum number of segments present in any iteration in the range selected'''
        n_particles = self.get_summary_table()['n_particles']
        return n_particles[first_iter-1:last_iter].max()
                
    def total_segs_in_range(self, first_iter, last_iter):
        '''Return the total number of segments present in all iterations in the range selected'''
        n_particles = self.get_summary_table()['n_particles']
        return n_particles[first_iter-1:last_iter].sum()

    def get_pcoord_len(self, n_iter):
        pcoord_ds = self.get_pcoord_dataset(n_iter)
        return pcoord_ds.shape[1]
    