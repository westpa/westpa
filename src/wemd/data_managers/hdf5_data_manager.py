from wemd.core import Segment, WESimIter
from wemd.data_managers import DataManagerBase
import h5py
from h5py.h5e import SymbolError
import numpy

from numpy import uint8, uint64, float64

import logging
log = logging.getLogger(__name__)

class SegIndexRow(object):
    dtype = numpy.dtype([('seg_id', uint64),
                         ('p_parent_id', uint64),
                         ('weight', float64),
                         ('cputime', float64),
                         ('walltime', float64),
                         ('status', uint8),])
    idx_by_field = dict((field, i) for (i,field) in enumerate(dtype.fields.iterkeys()))
    
    def __init__(self, item = None):
        self.data = numpy.zeros((1,), dtype=self.dtype)
        if item is not None:
            self.data[0] = item
    
    def __getitem__(self, key):
        return self.data[key]
    
    def __setitem__(self, key, value):
        self.data[key] = value
                
class HDF5DataManager(DataManagerBase):
    def __init__(self, runtime_config):
        super(HDF5DataManager,self).__init__(runtime_config)
        runtime_config.require('data.hdf5.filename')
        
        self.hdf5file = h5py.File(runtime_config['data.hdf5.filename'])
        try:
            self._iter_nodename_prec = self.hdf5file.attrs['iter_nodename_prec']
        except KeyError:
            self._iter_nodename_prec = runtime_config.get_int('data.hdf5.iter_nodename_prec', 12)
            self.hdf5file.attrs['nodename_prec'] = uint8(self._iter_nodename_prec)
        
        self.iterations_index_dtype = numpy.dtype([('nodename', 
                                                    numpy.dtype('|S%d' % self._iter_nodename_prec))])        
    def _iter_nodename_from_id(self, id_):
        return '%0.*d' % (self._iter_nodename_prec, id_)
            
    def _get_we_sim_iter_path(self, we_sim_iter):
        try:
            n_iter = we_sim_iter.n_iter
        except AttributeError:
            n_iter = we_sim_iter
            
        iteridx = self.hdf5file['/WE/Iterations/IterIndex']
        return '/WE/Iterations/%s' % (iteridx[n_iter, 'nodename'])    
                
    def prepare_backing(self, sim_config):
        self.hdf5file.create_group('/WE')
        iterations_group = self.hdf5file.create_group('/WE/Iterations')
        iterations_group.create_dataset('IterIndex', 
                                        shape = (1,),
                                        maxshape = (None,),
                                        chunks = (4,),
                                        dtype = self.iterations_index_dtype)
        
    def create_we_sim_iter(self, we_sim_iter):
        assert we_sim_iter.n_iter is not None
        assert we_sim_iter.n_particles > 0
        assert we_sim_iter.binarray is not None
        
        n_iter = we_sim_iter.n_iter
        n_particles = we_sim_iter.n_particles
        pcoord_ndim = we_sim_iter.binarray.ndim
                
        iterations_group = self.hdf5file['/WE/Iterations']
        idx = iterations_group['IterIndex']
        
        if n_iter > 0:
            # We are extending the simulation
            assert n_iter == len(idx)
            log.debug('extending index')
            idx.resize((n_iter+1,))
            
        iter_node_name = self._iter_nodename_from_id(n_iter)
        idx[n_iter] = (iter_node_name,)
        new_iter_group = iterations_group.create_group(iter_node_name)

        # Create the segment index
        segidx = new_iter_group.create_dataset('SegIndex', 
                                               shape = (n_particles,),
                                               dtype = SegIndexRow.dtype)
        
        # Create the progress coordinate array
        # TODO: add support for different data type
        new_iter_group.create_dataset('ProgCoords',
                                      shape=(n_particles, 1, pcoord_ndim),
                                      chunks=(2, 1, pcoord_ndim),
                                      maxshape=(n_particles, None, pcoord_ndim),
                                      dtype=float64)
        
        # Create the family tree array
        new_iter_group.create_dataset('SegmentLineage',
                                      shape=(n_particles, 1),
                                      dtype=uint64)
        
        self._update_we_sim_iter(we_sim_iter, new_iter_group)
        
    def _update_we_sim_iter(self, we_sim_iter, iter_group):
        if we_sim_iter.binarray is not None:
            bin_boundaries = numpy.array(we_sim_iter.binarray.boundaries)
            try:
                ds = iter_group.require_dataset('BinBoundaries', 
                                                shape=bin_boundaries.shape,
                                                dtype=bin_boundaries.dtype)
            except ValueError:
                log.warning('unexpected resize of BinBoundaries')
                del iter_group['BinBoundaries']
                ds = iter_group.create_dataset('BinBoundaries',
                                               shape=bin_boundaries.shape,
                                               dtype=bin_boundaries.dtype)
            ds[...] = bin_boundaries

            for (nodename, contents) in (('BinIdealNumParticles', we_sim_iter.binarray.ideal_num),
                                         ('BinSplitThreshold', we_sim_iter.binarray.split_threshold),
                                         ('BinMergeThresholdMin', we_sim_iter.binarray.merge_threshold_min),
                                         ('BinMergeThresholdMax', we_sim_iter.binarray.merge_threshold_max)):
                try:
                    ds = iter_group.require_dataset(nodename, shape=contents.shape, dtype=contents.dtype)
                except ValueError:
                    log.warning('unexpected resize of %s' % nodename)
                    del iter_group[nodename]
                    ds = iter_group.create_dataset(nodename, shape=contents.shape, dtype=contents.dtype)
                ds[...] = contents
            
            for (nodename, generating_func) in (('BinInitialPopulation', we_sim_iter.binarray.population_array),
                                                ('BinInitialNumParticles', we_sim_iter.binarray.nparticles_array)):
                contents = generating_func()
                try:
                    ds = iter_group.require_dataset(nodename, shape=contents.shape, dtype=contents.dtype)
                except ValueError:
                    log.warning('unexpected resize of %s' % nodename)
                    del iter_group[nodename]
                    ds = iter_group.create_dataset(nodename, shape=contents.shape, dtype=contents.dtype)
                ds[...] = contents
        else:
            for key in ('BinBoundaries', 'BinIdealNumParticles', 
                        'BinSplitThreshold', 'BinMergeThresholdMin',
                        'BinMergeThresholdMax', 'BinInitialPopulation',
                        'BinInitialNumParticles'):
                try:
                    del iter_group[key]
                except (KeyError,SymbolError):
                    pass

        iter_group.attrs['n_iter'] = uint64(we_sim_iter.n_iter)
        iter_group.attrs['norm'] = float64(we_sim_iter.norm)        
        iter_group.attrs['n_particles'] = uint64(we_sim_iter.n_particles)
        for (attr,dtype) in (('cputime', float64), ('walltime', float64)):
            val = getattr(we_sim_iter, attr)
            if val is None:
                try:
                    del iter_group.attrs[attr]
                except KeyError:
                    pass
            else:
                iter_group.attrs[attr] = dtype(val)

    def update_we_sim_iter(self, we_sim_iter):
        self._update_we_sim_iter(we_sim_iter, self.hdf5file[self._get_we_sim_iter_path(we_sim_iter)])

                
    def create_segments(self, we_sim_iter, segments):
        iter_group = self.hdf5file[self._get_we_sim_iter_path(we_sim_iter)]
        segidx = iter_group['SegIndex']
        pcoords = iter_group['ProgCoords']
        lineage = iter_group['SegmentLineage']
        
        # Check for duplicate IDs
        existing_seg_ids = set(segidx[:,'seg_id'])
        existing_seg_ids.discard(0)
        new_seg_ids = set(segment.seg_id for segment in segments)
        new_seg_ids.difference_update((None,0))
        duplicate_ids = existing_seg_ids & new_seg_ids
        if duplicate_ids:
            raise ValueError('duplicate key(s): %s' % list(duplicate_ids))
        
        # Find the maximum existing ID
        try:
            max_seg_id = max(existing_seg_ids)
        except ValueError:
            # Empty sequence
            max_seg_id = 0
        
        # Free some memory if we need to
        del existing_seg_ids, new_seg_ids
        
        # Find the first empty index row
        first_empty_row = 0
        for irow in xrange(0, len(segidx)):
            if segidx[irow,'seg_id'] == 0:
                first_empty_row = irow
                break
                
        for (inewseg, segment) in enumerate(segments):
            irow = first_empty_row + inewseg

            if segment.seg_id is None:
                seg_id = segment.seg_id = max_seg_id + inewseg + 1
            
            row = SegIndexRow(segidx[irow])
            row['seg_id'] = segment.seg_id
            segidx[irow] = row.data
        
        self.update_segments(we_sim_iter, segments)
                
            
    def update_segments(self, we_sim_iter, segments):
        iter_group = self.hdf5file[self._get_we_sim_iter_path(we_sim_iter)]
        segidx = iter_group['SegIndex']
        pcoords = iter_group['ProgCoords']
        lineage = iter_group['SegmentLineage']
        
        # Resize progress coordinate array if necessary
        shape_set = set((segment.pcoord is not None) and segment.pcoord.shape or None
                        for segment in segments)
        if len(shape_set) > 1:
            # Refuse to commit changes where segments have different
            # pcoord array sizes
            raise ValueError('pcoord shapes do not match')
        pcoord_shape = shape_set.pop()
        if pcoord_shape is not None and pcoords.shape[1:] != pcoord_shape:
            log.debug('resizing pcoord array to %r' % (pcoord_shape,))
            pcoords.resize(pcoords.shape[0] + pcoord_shape)

        # Resize the family tree if necessary
        max_n_parents = max(len(seg.parents) for seg in segments)
        if lineage.shape[1] < max_n_parents:
            log.debug('resizing particle lineage table')
            old_lineage = iter_group.copy(lineage, '_SegmentLineage')
            lineage = iter_group.create_dataset('SegmentLineage',
                                                shape=(old_lineage.shape[0],
                                                       max_n_parents),
                                                dtype=old_lineage.dtype)
            lineage[:, 0:old_lineage.shape[1]] = old_lineage
            del iter_group['_SegmentLineage']
            
        # Create an index into the index :)
        # i.e. have a mapping at hand between seg_id and array indices
        index_by_id = dict((segidx[idx, 'seg_id'], idx) for idx in xrange(0, len(segidx)))
        
        # Save updates
        for segment in segments:
            assert segment.n_iter == we_sim_iter.n_iter
            
            if segment.seg_id is None:
                raise ValueError('seg_id cannot be None')
            
            irow = index_by_id[segment.seg_id]
            idxrow = SegIndexRow(segidx[irow])            
            idxrow['weight'] = segment.weight
            idxrow['cputime'] = segment.cputime or 0.0
            idxrow['walltime'] = segment.walltime or 0.0
            idxrow['status'] = segment.status
            
            if segment.p_parent is not None:
                idxrow['p_parent_id'] = segment.p_parent.seg_id
            else:
                idxrow['p_parent_id'] = 0
                
            if segment.parents:
                lineage[irow] = [p.seg_id for p in segment.parents]
            else:
                # Broadcasts across entire row, setting all entries to zero
                lineage[irow] = 0
            
            if segment.pcoord is not None:
                pcoords[irow] = segment.pcoord
            else:
                pcoords[irow] = 0
                
            segidx[irow] = idxrow.data
                        
    def num_incomplete_segments(self, we_iter):
        ninc = 0
        segidx = self.hdf5file['%s/SegIndex' % self._get_we_sim_iter_path(we_iter)]
        for irow in xrange(0, len(segidx)):
            if segidx[irow, 'status'] != Segment.SEG_STATUS_COMPLETE:
                ninc += 1
        return ninc
    
    def get_we_sim_iter(self, n_iter):
        iter_node = self.hdf5file[self._get_we_sim_iter_path(n_iter)]
        we_iter = WESimIter(n_iter = iter_node.attrs['n_iter'],
                            norm   = iter_node.attrs['norm'],
                            n_particles = iter_node.attrs['n_particles'])
        assert we_iter.n_iter == n_iter
        return we_iter
    
    def get_segments(self, we_iter, seg_load_pcoord = True):
        iter_group = self.hdf5file[self._get_we_sim_iter_path(we_iter)]
        segments = []
        segidx = iter_group['SegIndex']
        for irow in xrange(0,len(segidx)):
            seg_id = segidx[irow,'seg_id']
            segment = Segment(seg_id = seg_id,
                              n_iter = iter_group.attrs['n_iter'],
                              status = segidx[irow, 'status'],
                              weight = segidx[irow, 'weight'],
                              cputime = segidx[irow, 'cputime'],
                              walltime = segidx[irow, 'walltime'],
                              )
            if seg_load_pcoord:
                segment.pcoord = iter_group['ProgCoords'][irow]
            segments.append(segment)
        return segments
    
    def get_prepared_segments(self, we_iter):
        iter_group = self.hdf5file[self._get_we_sim_iter_path(we_iter)]
        segments = []
        segidx = iter_group['SegIndex']
        for irow in xrange(0,len(segidx)):
            if segidx[irow, 'status'] != Segment.SEG_STATUS_PREPARED:
                continue
            seg_id = segidx[irow,'seg_id']
            segment = Segment(seg_id = seg_id,
                              n_iter = iter_group.attrs['n_iter'],
                              status = segidx[irow, 'status'],
                              weight = segidx[irow, 'weight'],
                              cputime = segidx[irow, 'cputime'],
                              walltime = segidx[irow, 'walltime'],
                              )
            segments.append(segment)
        return segments
        
        
