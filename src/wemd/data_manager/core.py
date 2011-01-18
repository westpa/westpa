from copy import copy
import logging
log = logging.getLogger(__name__)
import time
try:
    import cPickle as pickle
except ImportError:
    import pickle
    
import numpy
import h5py
import array

from wemd.core import Segment, WESimIter, Trajectory

def _n_iter(we_sim_iter):
    try:
        return we_sim_iter.n_iter
    except AttributeError:
        return we_sim_iter
    
class DataManagerBase(object):
    def __init__(self, runtime_config):
        self.runtime_config = runtime_config

class HDF5_data_manager(DataManagerBase):
    def __init__(self, runtime_config):
        super(HDF5_data_manager,self).__init__(runtime_config)
        runtime_config.require('data.db.h5')
        self.h5_name = runtime_config['data.db.h5']
        self.h5 = None
        self.index = None
        self.warn = None
        
    def prepare_database(self):
        self.h5 = h5py.File(self.h5_name)
    
    def require_hdf5(self):
        if self.h5 is None:
            self.h5 = h5py.File(self.h5_name)
        return self.h5
    
    def close_hdf5(self):
        if self.h5 is not None:
            self.h5.close()
            self.h5 = None
            
    def flush_hdf5(self):
        if self.h5 is not None:
            self.h5.flush()
            
    def create_we_sim_iter(self, we_sim_iter):
        h5 = self.require_hdf5()
        n_iter = _n_iter(we_sim_iter)

        cur_iter = h5.require_group('iter_%d' % n_iter)
        
        we_sim_tuple = ( we_sim_iter.n_iter,
                         we_sim_iter.n_particles,
                         we_sim_iter.norm,
                         we_sim_iter.binarray,
                         we_sim_iter.data)
        
        we_sim_data =  pickle.dumps(we_sim_tuple, pickle.HIGHEST_PROTOCOL)
        
        try:
            del cur_iter['we_sim_data']
        except KeyError:
            pass
        
        h5_we_sim_data = cur_iter.create_dataset('we_sim_data', data = we_sim_data)
        self.flush_hdf5()
        
    def update_we_sim_iter(self, we_sim_iter):
        self.create_we_sim_iter(we_sim_iter)
        
    def get_we_sim_iter(self, n_iter):  
        h5 = self.require_hdf5()        
        cur_iter = h5['iter_%d' % n_iter]
        we_sim_data = cur_iter['we_sim_data']       
        we_sim_data_str = we_sim_data[...]
        
        try:
            ( n_iter,
              n_particles,
              norm,
              binarray,
              data ) = pickle.loads(we_sim_data_str)
        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling data for WE iter %d' % n_iter)
            return None
        
        we_sim_iter =  WESimIter(n_iter = n_iter,
                         n_particles = n_particles,
                         norm = norm,
                         data = data)
        
        we_sim_iter.binarray = binarray
      
        return we_sim_iter
    
    def del_iter(self, n_iter):
        h5 = self.require_hdf5()    
            
        try:
            del h5['iter_%d' % n_iter]
        except:
            self.close_hdf5()
            raise
        self.flush_hdf5()
        
            
    def update_segments(self, we_sim_iter, segments):
       
        h5 = self.require_hdf5()
        n_iter = _n_iter(we_sim_iter)
        cur_iter = h5['iter_%d' % n_iter]
                
        sorted_segs = sorted(segments, key = (lambda s: s.seg_id) )
        
        min_seg_id = cur_iter['min_seg_id_nsegs'][0]
        nsegs = cur_iter['min_seg_id_nsegs'][1]
        
        sorted_seg_ids = [s.seg_id for s in sorted_segs]
        
        if None in sorted_seg_ids:
            raise ValueError('Cannot update segment without id')
        
        create_pcoord_db = False
        try:
            pcoord_shape = cur_iter['pcoord'].shape
        except KeyError:
            #if the original segments did not have pcoords, create pcoord now
            create_pcoord_db = True
        
        max_pcoord_len = 0
        ndim = 0
        #update pcoords
        for s in sorted_segs:
            if s.pcoord is not None:
                max_pcoord_len = len(s.pcoord)
                ndim = len(s.pcoord[0])
                            
        if max_pcoord_len != 0:
               
            if create_pcoord_db:
                h5_pcoord = cur_iter.create_dataset('pcoord', shape = (nsegs, max_pcoord_len, ndim), dtype = numpy.float64)               
            else:
                h5_pcoord = cur_iter['pcoord']
                
            pcoord_shape = h5_pcoord.shape
            #the code does not support variable length pcoords
            assert(max_pcoord_len == pcoord_shape[1])      
                
            i = 0
            for id in sorted_seg_ids:
                seg_off = id - min_seg_id
                assert(seg_off >= 0)      
                          
                if sorted_segs[i].pcoord is not None:
                    h5_pcoord[seg_off] = sorted_segs[i].pcoord
                i += 1
    
        #update weight
        h5_weight = cur_iter['weight']
        
        i = 0
        for id in sorted_seg_ids:
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
            
            h5_weight[seg_off] = sorted_segs[i].weight
            i += 1
            
        #update family tree
        max_parents_len = 0
        for i in xrange(0, len(sorted_segs)):
            n_parents = len(sorted_segs[i].parents)
            if n_parents > max_parents_len:
                max_parents_len = n_parents

        parents_shape = cur_iter['parents'].shape
        
        if max_parents_len + 3 > parents_shape[1]:
            cur_iter['parents'].resize((parents_shape[0], max_parents_len + 3))
        
        h5_parents = cur_iter['parents']
        
        #seg_id, p_parent, n_parents, parent ids
        for i in xrange(0, len(sorted_seg_ids)):
            id = sorted_seg_ids[i]
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
            
            h5_parents[seg_off,0] = sorted_seg_ids[i]
            
            if sorted_segs[i].p_parent is not None:
                p_parent_id = sorted_segs[i].p_parent.seg_id
            else:
                p_parent_id = 0
                
            h5_parents[seg_off,1] = p_parent_id
            
            if p_parent_id == 0 and len(sorted_segs[i].parents) != 0:
                if self.warn is None:
                    log.warning('p_parent is None, but parents are set -> Warning will not be repeated')
                    self.warn = 1
            #only update parents if not given during creation
            #segments sent to update segments usually do not have parents loaded
            parent_ids = [p.seg_id for p in sorted_segs[i].parents]
            
            if h5_parents[seg_off,2] == 0:
                len_parent_ids = len(parent_ids)
                if len_parent_ids > 0:
                    h5_parents[seg_off,2] = len_parent_ids
                    h5_parents[seg_off,3:len_parent_ids+3] = parent_ids
                        
        #seg_status
        i = 0
        for id in sorted_seg_ids:
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
            status = sorted_segs[i].status
            
            if not status:
                status = 0
            
            cur_iter['status'][seg_off] = status          
            i += 1
            
        #seg endpoint type
        i = 0
        for id in sorted_seg_ids:
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
            cur_iter['endpoint_type'][seg_off] = sorted_segs[i].endpoint_type            
            i += 1
        
        #timing info
        old_timing_pickle = list(cur_iter['timing'][...])

        i = 0
        for id in sorted_seg_ids:
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
     
            old_timing_pickle[seg_off] = pickle.dumps((sorted_segs[i].cputime, sorted_segs[i].walltime, sorted_segs[i].starttime, sorted_segs[i].endtime), pickle.HIGHEST_PROTOCOL)           
            i += 1

        del cur_iter['timing']
        cur_iter.create_dataset('timing', data=old_timing_pickle)
        
        #data/data_ref
        #need to copy the old data, add new data and rewrite table
        #otherwise strings are truncated
        old_data_pickle = list(cur_iter['seg_data'][...])

        i = 0
        for id in sorted_seg_ids:
            seg_off = id - min_seg_id
            assert(seg_off >= 0)
            
            old_data_pickle[seg_off] = pickle.dumps((sorted_segs[i].data, sorted_segs[i].data_ref), pickle.HIGHEST_PROTOCOL) 
            i += 1

        del cur_iter['seg_data']
        h5_seg_data = cur_iter.create_dataset('seg_data', data = old_data_pickle)
        
        self.flush_hdf5()
    
    def create_segments(self, we_sim_iter, segments, delete_old_segs = True):

        h5 = self.require_hdf5()
        n_iter = _n_iter(we_sim_iter)
        cur_iter = h5['iter_%d' % n_iter]
        
        if delete_old_segs:
            for key in ('pcoord', 'weight', 'parents', 'status', 'endpoint_type', 'timing', 'seg_data', 'min_seg_id_nsegs'):
                try:
                    del cur_iter[key]
                except KeyError:
                    pass
                        
        #sort segments by segid
        sorted_segs = sorted(segments, key = (lambda s: s.seg_id) )
        min_seg_id = sorted_segs[0].seg_id

        if n_iter >= 2:
            if min_seg_id is None:
                prev_iter = h5['iter_%d' % (n_iter - 1)]
                prev_min_seg_id_nsegs = prev_iter['min_seg_id_nsegs'][...]
                min_seg_id = prev_min_seg_id_nsegs[0] + prev_min_seg_id_nsegs[1]
        else:
            if min_seg_id is None:
                min_seg_id = 1
            
        cur_iter.create_dataset('min_seg_id_nsegs', data = numpy.array([min_seg_id, len(sorted_segs)], dtype=numpy.uint64), shape = (2,), dtype = numpy.uint64)
        
        pcoord_len = 0
        for seg in sorted_segs:
            if seg.pcoord is not None:
                pcoord_len = len(seg.pcoord)
                ndim = len(seg.pcoord[0])
                
        if pcoord_len != 0:
            h5_pcoord = cur_iter.create_dataset('pcoord', shape=(pcoord_len, ndim), dtype = numpy.float64)
            
            i = 0
            for seg in sorted_segs:
                if seg.pcoord is not None:
                    h5_pcoord[i] = seg.pcoord
                i += 1        
                
        #weight
        h5_weight = cur_iter.create_dataset('weight', shape=(len(sorted_segs),), dtype=numpy.float64)

        for i in xrange(0,len(sorted_segs)):
            h5_weight[i] = sorted_segs[i].weight
            
        #family tree
        max_parents_len = 0
        for i in xrange(0, len(sorted_segs)):
            n_parents = len(sorted_segs[i].parents)
            if n_parents > max_parents_len:
                max_parents_len = n_parents

        #seg_id, p_parent, n_parents, parent ids
        h5_parents = cur_iter.create_dataset('parents', shape=(len(sorted_segs), max_parents_len + 3), dtype=numpy.uint64, compression='gzip')
        
        for i in xrange(0, len(sorted_segs)):
            s = sorted_segs[i]
            #seg_id
            h5_parents[i,0] = s.seg_id or 0
            
            #p_parent
            if sorted_segs[i].p_parent is not None:
                p_parent_id = s.p_parent.seg_id
            else:
                p_parent_id = 0
                
            h5_parents[i,1] = p_parent_id
            
            #n_parents
            n_parents = len(s.parents)
            h5_parents[i,2] = n_parents
    
            if p_parent_id == 0 and n_parents != 0:
                if self.warn is None:
                    log.warning('p_parent is None, but parents are set -> Warning will not be repeated')
                    self.warn = 1
                    
            #parent_ids
            if n_parents > 0:
                h5_parents[i,3:n_parents+3] = [p.seg_id for p in s.parents] 
                                
        #seg_status
        h5_status = cur_iter.create_dataset('status', shape=(len(sorted_segs),), dtype=numpy.uint8)
        
        for i in xrange(0,len(sorted_segs)):
            if sorted_segs[i].status:
                h5_status[i] = sorted_segs[i].status
            else:
                h5_status[i] = 0
    
        #seg endpoint type
        h5_endpoint = cur_iter.create_dataset('endpoint_type', shape=(len(sorted_segs),), dtype=numpy.uint8)

        for i in xrange(0,len(sorted_segs)):
            h5_endpoint[i] = sorted_segs[i].endpoint_type
                    
        #timing info
        time_pickle = [pickle.dumps((s.cputime, s.walltime, s.starttime, s.endtime), pickle.HIGHEST_PROTOCOL) for s in sorted_segs]
        cur_iter.create_dataset('timing', data=time_pickle)        
        
        #data/data_ref
        data_pickle = [pickle.dumps((s.data, s.data_ref), pickle.HIGHEST_PROTOCOL) for s in sorted_segs]
        cur_iter.create_dataset('seg_data', data=data_pickle)

        self.flush_hdf5()
        
    def num_incomplete_segments(self, we_iter = None):
        h5 = self.require_hdf5()
        if not we_iter:
            we_iter = self.get_last_iter()
            
        n_iter = _n_iter(we_iter)
        cur_iter = h5['iter_%d' % n_iter]
        iter = self.get_we_sim_iter( n_iter )
        n_particles = iter.n_particles
        
        try:
            status = cur_iter['status']
        except KeyError:
            if n_iter == 0:
                return 0
            else:
                return n_particles
        
        num_segs = 0
        for i in xrange(0, len(status)):
            if status[i] != Segment.SEG_STATUS_COMPLETE:
                num_segs += 1

        return num_segs

    def get_segment_by_dataset(self, dataset, seg_id, seg_offset, load_p_parent = False):
        try:
            seg_pcoord = dataset['pcoord'][seg_offset][...]
        except KeyError:
            pcoord_len = 0
            seg_pcoord = None
            
        seg_weight = dataset['weight'][seg_offset]
        
        seg_p_parent_id = dataset['parents'][seg_offset,1]
        if seg_p_parent_id == 0:
            seg_p_parent_id = None
            
        '''
        if seg_p_parent_id is not None:
            len_parents = dataset['parents'][seg_offset,2]
            
            if len_parents > 0:
                seg_parents = dataset['parents'][seg_offset,3:len_parents+3]              
            else:
                seg_parents = None
        Currently we do not load parents, so there is no need to read parent info
        '''
        
        seg_status = dataset['status'][seg_offset]
        if seg_status == 0:
            seg_status = Segment.SEG_STATUS_UNSET
            
        seg_endpoint_type = dataset['endpoint_type'][seg_offset]
        
        seg_timing_pickle = dataset['timing'][seg_offset]

        try:
            seg_cputime, seg_walltime, seg_starttime, seg_endtime = pickle.loads(seg_timing_pickle)
        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling seg timing for seg %r' % seg_id)
            raise

        seg_data_pickle = dataset['seg_data'][seg_offset]
                
        try:
            seg_data, seg_data_ref = pickle.loads(seg_data_pickle)
        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling seg data for seg %r' % seg_id)
            raise
                    
        we_sim_data = dataset['we_sim_data']       
        we_sim_data_str = we_sim_data[...]
        
        try:
            ( n_iter,
              n_particles,
              norm,
              binarray,
              data ) = pickle.loads(we_sim_data_str)
        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling we data')
            raise
                          
        if load_p_parent and seg_p_parent_id:
            p_parent = self.get_segment(seg_p_parent_id, load_p_parent = False)
        else:
            p_parent = None
            
        return Segment( seg_id = seg_id,
                        n_iter = n_iter,
                        status = seg_status,
                        p_parent = p_parent,
                        endpoint_type = seg_endpoint_type,
                        weight = seg_weight,
                        pcoord = copy(seg_pcoord),
                        data_ref = seg_data_ref,  
                        starttime = seg_starttime,
                        endtime = seg_endtime,
                        cputime = seg_cputime,
                        walltime = seg_walltime,                              
                        data = seg_data
                        )
    
    def build_index(self):
        h5 = self.require_hdf5()
        
        self.index = []
        for iter in sorted(h5.keys()):
            i = int(iter.split('_')[1])
            
            if i == 0:
                continue
            
            min_seg_id_nsegs = h5[iter]['min_seg_id_nsegs'][:]
            
            self.index.insert(0, dict(iter = i,
                                   min_seg_id = min_seg_id_nsegs[0],
                                   nsegs = min_seg_id_nsegs[1] ))
            
    def update_index(self, max_iter):
        h5 = self.require_hdf5()
        
        if self.index is None:
            self.build_index()
            
        last_iter = self.index[0]['iter']
        for i in range(last_iter + 1, max_iter + 1):
            iter = ('iter_%d' % i)
            
            try:
                min_seg_id_nsegs = h5[iter]['min_seg_id_nsegs'][:]
            except KeyError:
                raise ValueError('invalid iter %d' % i)
            
            self.index.insert(0, dict(iter = i,
                                   min_seg_id = min_seg_id_nsegs[0],
                                   nsegs = min_seg_id_nsegs[1] ))    
                    
    def get_max_seg_id(self, we_iter):
        h5 = self.require_hdf5()
        n_iter = _n_iter(we_iter)
        if n_iter == 0:
            return 0
        
        cur_iter = h5['iter_%d' % n_iter]
        ids = cur_iter['min_seg_id_nsegs'][:]
        return int(ids[0] + ids[1] - 1)
    
    def get_last_iter(self):
        h5 = self.require_hdf5()
        
        iters = [int(k.split('_')[1]) for k in h5.keys()]
        sorted_iter = sorted(iters)
        max_iter = sorted_iter[-1]       

        return self.get_we_sim_iter(max_iter)

    def get_last_complete_iter(self):
        h5 = self.require_hdf5()
        
        iters = [int(k.split('_')[1]) for k in h5.keys()]
        sorted_iter = sorted(iters)
        max_iter = sorted_iter[-1]       

        return self.get_we_sim_iter(max_iter-1)
    
    def get_first_iter(self):
        h5 = self.require_hdf5()
        
        iters = [int(k.split('_')[1]) for k in h5.keys()]
        sorted_iter = sorted(iters)
        max_iter = max(sorted_iter[0],1) #iteration 0 does not contain any segments       

        return self.get_we_sim_iter(max_iter)    
        
    def get_segment(self, seg_id, load_p_parent = False):
        h5 = self.require_hdf5()

        iters = [int(k.split('_')[1]) for k in h5.keys()]
        sorted_iter = sorted(iters)
        max_iter = sorted_iter[-1]
               
        if self.index is None:
            self.build_index()
        elif self.index[-1]['iter'] < max_iter:
            self.update_index(max_iter)

        for i in self.index:
            min_seg_id = i['min_seg_id']
            nsegs = i['nsegs']
            
            if seg_id >= min_seg_id and seg_id < (min_seg_id + nsegs):
                seg_offset = seg_id - min_seg_id
                iter = 'iter_%d' % i['iter']
                return self.get_segment_by_dataset(h5[iter], seg_id, seg_offset, load_p_parent = load_p_parent)
            
        log.warning('could not find segment seg_id %r' % seg_id)
        return None
    
    def get_segments(self, we_n_iter, n_iter_upper = None, status_criteria = None, status_negate = False, endpoint_criteria = None, endpoint_negate = False, load_p_parent = False):
        h5 = self.require_hdf5()
        
        if n_iter_upper is None:
            n_iter_upper = we_n_iter
            
        segs = []
        for iter_num in xrange(we_n_iter, n_iter_upper + 1):
            iter = h5['iter_%d' % iter_num]
            min_seg_id = iter['min_seg_id_nsegs'][0]
            nsegs = iter['min_seg_id_nsegs'][1]
            
            seg_offsets = range(0, nsegs)
            
            try:
                pcoords = iter['pcoord'][seg_offsets][:]
            except KeyError:
                pcoord_len = 0
                pcoords = None
                
            weights = iter['weight'][seg_offsets][:]
            
            p_parent_ids = iter['parents'][seg_offsets,1][:]
       
            '''
            for ...
            seg_parents = []
            if p_parent_id is not None:
                len_parents = iter['parents'][seg_offsets,2]
                for j in xrange(0, len(len_parents)):
                    if len_parents[j] > 0:
                        seg_parents.append(iter['parents'][seg_offsets[j],3:len_parents[j]+3])              
                    else:
                        seg_parents.append(None)
            '''
            
            status = iter['status'][seg_offsets][:]
                
            endpoint_types = iter['endpoint_type'][seg_offsets][:]
            
            timing_pickles = iter['timing'][seg_offsets][:]
            timing_values = [] # [(seg_cputime, seg_walltime, seg_starttime, seg_endtime),...]
            for timing_pickle in timing_pickles:
                try:
                    timing_values.append( pickle.loads(timing_pickle) )
                except (pickle.UnpicklingError, EOFError):
                    log.error('Error unpickling seg timing for seg %r' % seg_offsets)
                    raise
    
            data_pickles = iter['seg_data'][seg_offsets][:]
            data_values = [] #[(seg_data, seg_data_ref), ...]
            for data_pickle in data_pickles:                     
                try:
                    data_values.append(pickle.loads(data_pickle))
                except (pickle.UnpicklingError, EOFError):
                    log.error('Error unpickling seg data for seg %r' % seg_offsets)
                    raise
                        
            we_sim_data = iter['we_sim_data']       
            we_sim_data_str = we_sim_data[...]
            
            try:
                ( n_iter,
                  n_particles,
                  norm,
                  binarray,
                  data ) = pickle.loads(we_sim_data_str)
            except (pickle.UnpicklingError, EOFError):
                log.error('Error unpickling we data')
                raise
            
            p_parents = []
            
            if load_p_parent:
                for k in xrange(0,len(p_parent_ids)):
                    if p_parent_ids[k]:
                        p_parents.append(self.get_segment(p_parent_ids[k], load_p_parent = False))
                    else:
                        p_parents.append(None)
            else:
                p_parents = [None for k in p_parent_ids]
            
            
            for i in xrange(0, len(seg_offsets)): 
                               
                if status[i] == 0:
                    stat = Segment.SEG_STATUS_UNSET
                else:
                    stat = status[i]
                                            
                if status_criteria is not None:
                    if status_negate == False:
                        if status[i] != status_criteria:
                            continue
                    else:
                        if status[i] == status_criteria:
                            continue    

                if endpoint_criteria is not None:               
                    if endpoint_negate == False:
                        if endpoint_types[i] != endpoint_criteria:
                            continue
                    else:
                        if endpoint_types[i] == endpoint_criteria:
                            continue                                    
                         
                seg_data, seg_data_ref = data_values[i]
                seg_cputime, seg_walltime, seg_starttime, seg_endtime = timing_values[i]
                
                pcoord = None
                if pcoords is not None:
                    if pcoords[i] is not None:
                        pcoord = pcoords[i]
                        
                seg = Segment( seg_id = seg_offsets[i] + min_seg_id,
                            n_iter = n_iter,
                            status = stat,
                            p_parent = p_parents[i],
                            endpoint_type = endpoint_types[i],
                            weight = weights[i],
                            pcoord = copy(pcoord),
                            data_ref = seg_data_ref,  
                            starttime = seg_starttime,
                            endtime = seg_endtime,
                            cputime = seg_cputime,
                            walltime = seg_walltime,                              
                            data = seg_data
                            )

                segs.append(seg)
   
        return segs 

    def get_segments_by_parent_id(self, we_n_iter, n_iter_upper = None, p_parent_id = None, p_parent_id_negate = False):
        h5 = self.require_hdf5()
        
        if n_iter_upper is None:
            n_iter_upper = we_n_iter
            
        segs = []
        for iter_num in xrange(we_n_iter, n_iter_upper + 1):
            iter = h5['iter_%d' % iter_num]
            min_seg_id = iter['min_seg_id_nsegs'][0]
            nsegs = iter['min_seg_id_nsegs'][1]
            
            seg_offsets = range(0, nsegs)
                                    
            p_parent_ids = iter['parents'][seg_offsets,1][:]
            
            load_segs = []       
            for i in xrange(0,len(p_parent_ids)):
                if p_parent_id is None:
                    p_parent_id = 0
                
                if p_parent_id_negate:
                    if p_parent_ids[i] != p_parent_id:
                        load_segs.append(i)                    
                else:
                    if p_parent_ids[i] == p_parent_id:
                        load_segs.append(i)
                       
            if not load_segs:
                continue

            seg_offsets = load_segs
                         
            try:
                pcoords = iter['pcoord'][seg_offsets][:]
            except KeyError:
                pcoord_len = 0
                pcoords = None            
            
            weights = iter['weight'][seg_offsets][:]            
            status = iter['status'][seg_offsets][:]
                
            endpoint_types = iter['endpoint_type'][seg_offsets][:]
            
            timing_pickles = iter['timing'][seg_offsets][:]
            timing_values = [] # [(seg_cputime, seg_walltime, seg_starttime, seg_endtime),...]
            for timing_pickle in timing_pickles:
                try:
                    timing_values.append( pickle.loads(timing_pickle) )
                except (pickle.UnpicklingError, EOFError):
                    log.error('Error unpickling seg timing for seg %r' % seg_offsets)
                    raise
    
            data_pickles = iter['seg_data'][seg_offsets][:]
            data_values = [] #[(seg_data, seg_data_ref), ...]
            for data_pickle in data_pickles:                     
                try:
                    data_values.append(pickle.loads(data_pickle))
                except (pickle.UnpicklingError, EOFError):
                    log.error('Error unpickling seg data for seg %r' % seg_offsets)
                    raise
                        
            we_sim_data = iter['we_sim_data']       
            we_sim_data_str = we_sim_data[...]
            
            try:
                ( n_iter,
                  n_particles,
                  norm,
                  binarray,
                  data ) = pickle.loads(we_sim_data_str)
            except (pickle.UnpicklingError, EOFError):
                log.error('Error unpickling we data')
                raise
            
            p_parents = []        
            for k in seg_offsets:
                if p_parent_ids[k]:
                    p_parents.append(self.get_segment(p_parent_ids[k], load_p_parent = False))
                else:
                    p_parents.append(None)
            
            
            for i in xrange(0, len(seg_offsets)): 
                               
                if status[i] == 0:
                    stat = Segment.SEG_STATUS_UNSET
                else:
                    stat = status[i]                              
                         
                seg_data, seg_data_ref = data_values[i]
                seg_cputime, seg_walltime, seg_starttime, seg_endtime = timing_values[i]
                
                pcoord = None
                if pcoords is not None:
                    if pcoords[i] is not None:
                        pcoord = pcoords[i]
                        
                seg = Segment( seg_id = seg_offsets[i] + min_seg_id,
                            n_iter = n_iter,
                            status = stat,
                            p_parent = p_parents[i],
                            endpoint_type = endpoint_types[i],
                            weight = weights[i],
                            pcoord = copy(pcoord),
                            data_ref = seg_data_ref,  
                            starttime = seg_starttime,
                            endtime = seg_endtime,
                            cputime = seg_cputime,
                            walltime = seg_walltime,                              
                            data = seg_data
                            )

                segs.append(seg)
   
        return segs 

    def get_connectivity(self, we_n_iter, n_iter_upper = None):
        h5 = self.require_hdf5()

        segs = self.get_segments_by_parent_id(we_n_iter, n_iter_upper = n_iter_upper, p_parent_id = None, p_parent_id_negate = True)
        conn = [[s.p_parent.seg_id,s.seg_id] for s in segs]
        return conn