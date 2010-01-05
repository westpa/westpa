from wemd.core import Segment
from wemd.data_managers import DataManagerBase
import h5py
from h5py.h5e import SymbolError
import numpy

from numpy import uint8, uint64, float64

import logging
log = logging.getLogger(__name__)

ITERIDX_HDF5REF = 0
SEGIDX_HDF5REF = 0

class HDF5DataManager(DataManagerBase):
    def __init__(self, runtime_config):
        super(HDF5DataManager,self).__init__(runtime_config)
        runtime_config.require('data.hdf5.filename')
        
        self.hdf5file = h5py.File(runtime_config['data.hdf5.filename'])
        try:
            self._nodename_prec = self.hdf5file.attrs['nodename_prec']
        except KeyError:
            self._nodename_prec = runtime_config.get_int('data.hdf5.nodename_prec', 12)
            self.hdf5file.attrs['nodename_prec'] = uint8(self._nodename_prec)
        
        self.iterations_index_dtype = numpy.dtype([('hdf5ref', 
                                                    numpy.dtype('|S%d' % self._nodename_prec))])
        self.segments_index_dtype = numpy.dtype([('hdf5ref',
                                                  numpy.dtype('|S%d' % self._nodename_prec))])

        
    def _nodename_from_id(self, id_):
        return '%0.*d' % (self._nodename_prec, id_)
    
    def _get_we_sim_iter_path(self, we_sim_iter):
        try:
            we_iter = we_sim_iter.we_iter
        except AttributeError:
            we_iter = we_sim_iter
            
        iteridx = self.hdf5file['/Iterations/Index']
        return '/Iterations/%s' % (iteridx[we_iter][ITERIDX_HDF5REF])    
        
    def _get_segment_path(self, we_iter, seg_id):
        segpath = self._get_we_sim_iter_path(we_iter)
        segidx = self.hdf5file['%s/Segments/Index' % segpath]
        return '%s/Segments/%s' % (segpath, segidx[seg_id][SEGIDX_HDF5REF])
        
    def prepare_backing(self, sim_config):
        iterations_group = self.hdf5file.create_group('Iterations')
        iterations_group.create_dataset('Index', 
                                        shape = (1,),
                                        maxshape = (None,),
                                        #chunks = (1024,),
                                        dtype = self.iterations_index_dtype)
        
    def create_we_sim_iter(self, we_sim_iter):
        assert we_sim_iter.we_iter is not None
                
        iterations_group = self.hdf5file['Iterations']
        idx = iterations_group['Index']
        
        if we_sim_iter.we_iter == 0:
            # We are bootstrapping the simulation
            log.debug('bootstrapping simulation; not resizing index')
        else:
            # We are extending the simulation
            assert we_sim_iter.we_iter == len(idx)
            log.debug('extending index')
            idx.resize((len(idx)+1,))
        log.debug('index shape: %r' % idx.shape)
            
        iter_node_name = self._nodename_from_id(we_sim_iter.we_iter)
        idx[we_sim_iter.we_iter] = (iter_node_name,)
        new_iter_group = iterations_group.create_group(iter_node_name)

        seg_group = new_iter_group.create_group('Segments')
        segidx = seg_group.create_dataset('Index', 
                                        shape = (1,),
                                        maxshape = (None,),
                                        #chunks = (1024,),
                                        dtype = self.segments_index_dtype)
        self._update_we_sim_iter(we_sim_iter, new_iter_group)
        
    def _update_we_sim_iter(self, we_sim_iter, iter_group):
        if we_sim_iter.binarray is not None:
            try:
                iter_group.create_dataset('BinBoundaries', data=numpy.array(we_sim_iter.binarray.boundaries))
            except ValueError:
                iter_group['BinBoundaries'][...] = numpy.array(we_sim_iter.binarray.boundaries)
            for (nodename, contents) in (('BinIdealNumParticles', we_sim_iter.binarray.ideal_num),
                                         ('BinSplitThresholds', we_sim_iter.binarray.split_threshold),
                                         ('BinMergeThresholdMin', we_sim_iter.binarray.merge_threshold_min),
                                         ('BinMergeThresholdMax', we_sim_iter.binarray.merge_threshold_max)):
                log.debug('assigning to %r: %r' % (nodename, contents))
                try:
                    log.debug('shape of %r: %r' % (nodename, contents.shape))
                except AttributeError:
                    pass
                
                try:
                    iter_group.create_dataset(nodename, data=contents)
                except ValueError:
                    iter_group[nodename] = contents
            
            for (nodename, generating_func) in (('BinInitialPopulation', we_sim_iter.binarray.population_array),
                                                ('BinInitialNumParticles', we_sim_iter.binarray.nparticles_array)):
                contents = generating_func()
                try:
                    iter_group.create_dataset(nodename, data=contents)
                except ValueError:
                    iter_group[nodename] = contents

        iter_group.attrs['we_iter'] = uint64(we_sim_iter.we_iter)
        iter_group.attrs['norm'] = float64(we_sim_iter.norm)        
        iter_group.attrs['n_particles'] = uint64(we_sim_iter.n_particles or 0)
        iter_group.attrs['cputime'] = float64(we_sim_iter.cputime or 0.0)
        iter_group.attrs['walltime'] = float64(we_sim_iter.walltime or 0.0)
        
    
    def update_we_sim_iter(self, we_sim_iter):
        self._update_we_sim_iter(we_sim_iter, self.hdf5file[self._get_we_sim_iter_path(we_sim_iter)])
    
    def create_segment(self, segment):
        assert segment.seg_id is None
        
        segs_group = self.hdf5file['%s/Segments' % self._get_we_sim_iter_path(segment.we_iter)]
        segidx = segs_group['Index']
        last_seg_id = segidx.shape[0]-1
        if last_seg_id == 0:
            if segidx[0][0] == '':
                last_seg_id = None
        
        if last_seg_id is None:
            segment.seg_id = 0
        else:
            segidx.resize((len(segidx)+1,))
            segment.seg_id = last_seg_id+1
        segidx[segment.seg_id] = (self._nodename_from_id(segment.seg_id),)
            
        seg_group = segs_group.create_group(segidx[segment.seg_id][SEGIDX_HDF5REF])
        self._update_segment(segment, seg_group)
        
    def _update_segment(self, segment, seg_group):
        seg_group.attrs['seg_id'] = uint64(segment.seg_id)
        seg_group.attrs['we_iter'] = uint64(segment.we_iter)
        seg_group.attrs['weight'] = float64(segment.weight)
        seg_group.attrs['status'] = uint8(segment.status or 0)
        seg_group.attrs['cputime'] = float64(segment.cputime or 0.0)
        seg_group.attrs['walltime'] = float64(segment.walltime or 0.0)
        if segment.data_ref:
            seg_group.attrs['data_ref'] = numpy.array(segment.data_ref)
        else:
            try:
                del seg_group.attrs['data_ref']
            except KeyError:
                pass
        
        if segment.pcoord is not None:
            if len(segment.pcoord.shape) == 1:
                # Vector of a single n-dimensional progress coordinate
                pcoord_node = seg_group.require_dataset('Pcoord', dtype=segment.pcoord.dtype,
                                                       shape=(1, segment.pcoord.shape[0]),
                                                       )
                pcoord_node[0,:] = segment.pcoord
            else:
                # Array of time-resolved progress coordinates                
                try:
                    pcoord_node = seg_group.require_dataset('Pcoord', dtype=segment.pcoord.dtype,
                                                           shape = segment.pcoord.shape,
                                                           )
                except ValueError:
                    # Needs resize
                    del seg_group['Pcoord']
                    pcoord_node = seg_group.create_dataset('Pcoord', dtype=segment.pcoord.dtype,
                                                           shape = segment.pcoord.shape,
                                                           )
                pcoord_node[...] = segment.pcoord


        if segment.p_parent:
            seg_group.attrs['p_parent_id'] = uint64(segment.p_parent.seg_id)
        else:
            try:
                del seg_group.attrs['p_parent_id']
            except KeyError:
                pass
            
        if segment.parents:
            try:
                parents_node = seg_group.require_dataset('Parents', dtype=uint64,
                                                         shape = (len(segment.parents),))
                parents_node[...] = [p.seg_id for p in segment.parents]
            except ValueError:
                del seg_group['Parents']
                parents_node = seg_group.create_dataset('Parents', dtype=uint64,
                                                         shape = (len(segment.parents),))
                parents_node[...] = [p.seg_id for p in segment.parents]

    def update_segment(self, segment):
        self._update_segment(segment, self.hdf5file[self._get_segment_path(segment.we_iter, segment.seg_id)])
    
    def num_incomplete_segments(self, we_iter):
        ninc = 0
        seg_group = self.hdf5file['%s/Segments' % self._get_we_sim_iter_path(we_iter)]
        for row in seg_group['Index']:
            nodename = row[SEGIDX_HDF5REF]
            if seg_group[nodename].attrs['status'] != Segment.SEG_STATUS_COMPLETE:
                ninc += 1
        return ninc
    
    def get_segment(self, we_sim_iter, seg_id, load_p_parent = False):
        segment = Segment()
        seg_group = self.hdf5file[self._get_segment_path(we_sim_iter, seg_id)]
        segment.we_iter = seg_group.attrs['we_iter']
        segment.seg_id = seg_group.attrs['seg_id']
        segment.status = seg_group.attrs['status']
        segment.weight = seg_group.attrs['weight']
        segment.cputime = seg_group.attrs.get('cputime', 0.0)
        segment.walltime = seg_group.attrs.get('walltime', 0.0)
        segment.data_ref = seg_group.attrs.get('data_ref',None)
        try:
            segment.pcoord = seg_group['Pcoord'][...]
        except KeyError:
            segment.pcoord = None
        return segment
            
    def get_segments(self, we_sim_iter):
        sim_iter_path = self._get_we_sim_iter_path(we_sim_iter)
        segs_group = self.hdf5file['%s/Segments' % sim_iter_path]
        segments = [self.get_segment(we_sim_iter, seg_id)
                    for seg_id in xrange(0, segs_group['Index'].shape[0])]
        return segments
