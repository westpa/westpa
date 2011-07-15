from __future__ import print_function, division; __metaclass__=type
from itertools import imap, izip
import numpy
from wemd import Segment

import logging
log = logging.getLogger(__name__)

class CachingDataReader:
    def __init__(self, data_manager, cache_pcoords = False):
        self._data_manager = data_manager
        self._data_manager.open_backing(mode='r')
                
        # Mapping of iteration -> data from that iteration
        self._iter_groups = dict()
        self._seg_id_ranges = dict()
        self._seg_indices = dict()
        self._parent_arrays = dict()
        self._p_parent_arrays = dict()
        self._pcoord_arrays = dict()
        self._pcoord_datasets = dict()

        # Whether pcoord caching is being used
        self._cache_pcoords = cache_pcoords
                
        # For functions returning lists of segments, whether to include pcoords.
        self.include_pcoords = False
        
    def clear_cache(self):
        del self._iter_groups, self._seg_id_ranges, self._seg_indices, self._parent_arrays, self._p_parent_arrays
        del self._pcoord_arrays, self._pcoord_datasets
        
        self._iter_groups = dict()
        self._seg_id_ranges = dict()
        self._seg_indices = dict()
        self._parent_arrays = dict()
        self._p_parent_arrays = dict()
        self._pcoord_arrays = dict()
        self._pcoord_datasets = dict()
        
    
    def __getattr__(self, attr):
        return getattr(self._data_manager,attr)

    @property
    def cache_pcoords(self):
        '''Whether or not to cache progress coordinate data. While caching this data
        can significantly speed up some analysis operations, this requires
        copious RAM.
        
        Setting this to False when it was formerly True will release any cached data. 
        '''
        return self._cache_pcoords

    @cache_pcoords.setter
    def cache_pcoords(self, cache):
        self._cache_pcoords = cache
        
        if not cache:
            del self._pcoord_arrays
            self._pcoord_arrays = dict()
            
    def get_iter_group(self, n_iter):
        '''Return the HDF5 group corresponding to ``n_iter``'''
        try:
            return self._iter_groups[n_iter]
        except KeyError:
            iter_group = self._data_manager.get_iter_group(n_iter)
            return iter_group
    
    def get_segments(self, n_iter, include_pcoords = None):
        '''Return all segments present in iteration n_iter'''
        return self.get_segments_by_id(n_iter, self.get_seg_ids_where(n_iter, None), include_pcoords)
    
    def get_segments_by_id(self, n_iter, seg_ids, include_pcoords = None):
        '''Get segments from the data manager, employing caching where possible'''
        
        if include_pcoords is None: include_pcoords = self.include_pcoords
        
        if len(seg_ids) == 0: return []
        
        seg_index  = self.get_seg_index(n_iter)
        all_parent_ids = self.get_parent_array(n_iter)
        
        segments = []

        if include_pcoords:
            pcoords = self.get_pcoords(n_iter, seg_ids)
            #for (segment, pcoord_set) in izip(segments, pcoords):
            #    segment.pcoord = pcoord_set
        
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
    
    def get_children(self, segment):
        p_parents = self.get_p_parent_array(segment.n_iter+1)
        seg_ids = self.get_seg_ids_where(segment.n_iter+1, p_parents == segment.seg_id)
        return self.get_segments_by_id(segment.n_iter+1, seg_ids)

    def get_seg_index(self, n_iter):
        try:
            return self._seg_indices[n_iter]
        except KeyError:
            seg_index = self._seg_indices[n_iter] = self.get_iter_group(n_iter)['seg_index'][...]
            return seg_index
        
    def get_parent_array(self, n_iter):
        try:
            return self._parent_arrays[n_iter]
        except KeyError:
            parent_array = self._parent_arrays[n_iter] = self.get_iter_group(n_iter)['parents'][...]
            return parent_array
                
    def get_p_parent_array(self, n_iter):
        try:
            return self._p_parent_arrays[n_iter]
        except KeyError:
            parent_offsets = self.get_seg_index(n_iter)['parents_offset'][:]
            p_parent_array = self.get_parent_array(n_iter)[parent_offsets]
            self._p_parent_arrays[n_iter] = p_parent_array
            return p_parent_array
                
    def get_pcoord_array(self, n_iter):
        try:
            return self._pcoord_arrays[n_iter]
        except KeyError:
            pcoords = self._pcoord_arrays[n_iter] = self.get_iter_group(n_iter)['pcoord'][...]
            return pcoords
        
    def get_pcoord_dataset(self, n_iter):
        try:
            return self._pcoord_datasets[n_iter]
        except KeyError:
            pcoord_ds = self._pcoord_datasets[n_iter] = self.get_iter_group(n_iter)['pcoord']
            return pcoord_ds
        
    def get_pcoords(self, n_iter, seg_ids):
        if self._cache_pcoords:
            pcarray = self.get_pcoord_array(n_iter)
            return [pcarray[seg_id,...] for seg_id in seg_ids]
        else:
            return self.get_pcoord_dataset(n_iter)[list(seg_ids),...]
        
    def get_seg_ids_where(self, n_iter, bool_array = None):
        try:
            all_ids = self._seg_id_ranges[n_iter]
        except KeyError:
            all_ids = self._seg_id_ranges[n_iter] = numpy.arange(0,len(self.get_seg_index(n_iter)), dtype=numpy.uint32)
            
            
        if bool_array is None:             
            return all_ids
        else:        
            seg_ids = all_ids[bool_array]
            #seg_ids = [seg_id for (isegid, seg_id) in enumerate(all_ids) if bool_array[isegid]]        
            try:
                if len(seg_ids) == 0: return []
            except TypeError:
                # Not iterable
                return [seg_ids]
            else:
                return seg_ids
