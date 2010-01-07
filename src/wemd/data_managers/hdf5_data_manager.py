from wemd.core import Segment, WESimIter
from wemd.data_managers import DataManagerBase
import h5py
from h5py.h5e import SymbolError
import numpy

from numpy import uint8, uint64, float64

import logging
log = logging.getLogger(__name__)

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
        self.seg_index_dtype = numpy.dtype([('seg_id', uint64),
                                            ('p_parent_id', uint64),
                                            ('weight', float64),
                                            ('cputime', float64),
                                            ('walltime', float64),
                                            ('status', uint8),])
        
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
        
        if n_iter == 0:
            # We are bootstrapping the simulation
            log.debug('bootstrapping simulation; not resizing index')
        else:
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
                                               dtype = self.seg_index_dtype)
        
        # Create the progress coordinate array
        # TODO: add support for different data type
        new_iter_group.create_dataset('ProgCoords',
                                      shape=(n_particles, 1, pcoord_ndim),
                                      dtype=float64)
        
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
                
        self._update_segments(we_sim_iter, iter_group)
                
    def _update_segments(self, we_sim_iter, iter_group):
        segments = we_sim_iter.segments
        assert segments is not None and len(segments) > 0
        segidx = iter_group['SegIndex']
        idxrow = numpy.empty((1,), dtype=self.seg_index_dtype)
        
        pcoords = None
        if segments[0].pcoord is not None:
            pcoords = iter_group['ProgCoords']
            if pcoords[0].shape != segments[0].pcoord.shape:
                log.info('resizing progress coordinate storage')
                del iter_group['ProgCoords']
                pcoords = iter_group.create_dataset('ProgCoords',
                                                    shape = (we_sim_iter.n_particles,) + segments[0].pcoord.shape,
                                                    dtype = segments[0].pcoord.dtype,
                                                    chunks = (1,) + segments[0].pcoord.shape)
        
        segment_lineage = numpy.zeros((we_sim_iter.n_particles,1), dtype=uint64)
        for (iseg, segment) in enumerate(segments):
            # update index   
            idxrow[:] = segidx[iseg]
            if segment.seg_id is None:
                segment.seg_id = iseg + 1
                idxrow['seg_id'] = segment.seg_id
            else:
                idxrow['seg_id'] = segment.seg_id
            try:
                idxrow['p_parent_id'] = segment.p_parent.seg_id
            except AttributeError:
                idxrow['p_parent_id'] = 0
            idxrow['weight'] = segment.weight
            idxrow['cputime'] = segment.cputime or 0.0
            idxrow['walltime'] = segment.cputime or 0.0
            idxrow['status'] = segment.status
            segidx[iseg] = idxrow
            
            if pcoords is not None: 
                pcoords[iseg] = segment.pcoord
                
            if len(segment.parents) > segment_lineage.shape[1]:
                segment_lineage.resize((segment_lineage.shape[0], len(segment.parents)))
                segment_lineage[:,len(segment.parents)-1] = 0
            segment_lineage[iseg,0:len(segment.parents)] = [s.seg_id for s in segment.parents]
            
        # Store family tree
        try:
            del iter_group['SegmentLineage']
        except (KeyError,SymbolError):
            pass
        if (segment_lineage != 0).any():
            iter_group.create_dataset('SegmentLineage', data=segment_lineage)
            
                
        # Sanity checks
        assert set(segidx[:,'seg_id']) == set(xrange(1,len(segidx)+1))
            
    def update_we_sim_iter(self, we_sim_iter):
        self._update_we_sim_iter(we_sim_iter, self.hdf5file[self._get_we_sim_iter_path(we_sim_iter)])
        
    def num_incomplete_segments(self, we_iter):
        ninc = 0
        segidx = self.hdf5file['%s/SegIndex' % self._get_we_sim_iter_path(we_iter)]
        for irow in xrange(0, len(segidx)):
            if segidx[irow, 'status'] != Segment.SEG_STATUS_COMPLETE:
                ninc += 1
        return ninc
    
    def get_we_sim_iter(self, n_iter, load_segments = False):
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
        
        
