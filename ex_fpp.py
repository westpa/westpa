from __future__ import print_function, division; __metaclass__ = type
import sys, os, logging
import h5py, numpy
import wemd
from wt2.tool_classes import WEMDTool, WEMDDataReader, SegSelector
from wt2.tool_classes.data_reader import ByIterDataSelection
from wt2.tool_classes.selected_segs import SegmentSelection
from wemd.data_manager import n_iter_dtype, seg_id_dtype, weight_dtype
from trajtree import TrajTreeSet
import numpy, h5py

class FPPWalker:
    def __init__(self, h5file, aux_h5file):
        self.datasets = [ByIterDataSelection(aux_h5file, 'dssp', numpy.index_exp[:,87,'accessibility'], 'F19_sas', 'index'),
                         ByIterDataSelection(aux_h5file, 'dssp', numpy.index_exp[:,91,'accessibility'], 'W23_sas', 'index'),
                         ByIterDataSelection(aux_h5file, 'dssp', numpy.index_exp[:,94,'accessibility'], 'L26_sas', 'index'),
                         ]
        
        self.thresholds = [(1-0.946)*222.8, (1-0.876)*266.3, (1-0.865)*193.1]
        self.first_transitions = numpy.zeros((3,), n_iter_dtype)
        self.transition_weights = numpy.zeros((3,), weight_dtype)
        self.p_first = numpy.zeros((3,), weight_dtype)
        self.n_first = numpy.zeros((3,), numpy.uint)
                
    def visit(self, n_iter, seg_id, weight, has_children):
        sas_traces = numpy.array([dataset[n_iter, seg_id][:] for dataset in self.datasets], dtype=numpy.float32)
        ttimes = numpy.zeros((3,), numpy.int16)
        for (i, threshold) in enumerate(self.thresholds):
            under_threshold = sas_traces[i] < threshold
            idx_under = numpy.atleast_1d(numpy.squeeze(numpy.argwhere(under_threshold)))
            if idx_under.size:
                idx_under.sort()
                ttimes[i] = idx_under[0]
        if (ttimes > 0).any():
            first = ttimes.argmin()
            self.p_first[first] += weight
            self.n_first[first] += 1
            raise StopIteration
                    
data_manager = wemd.rc.get_data_manager()
data_manager.we_h5filename = 'system.h5'
data_manager.open_backing(mode='r')

aux_h5file = h5py.File('msms_dssp.h5', 'r')

fppw = FPPWalker(data_manager.we_h5file, aux_h5file)
segsel = SegmentSelection.from_text('binding_trajs.txt')
print('There are {:d} segments spanning {:d} iterations in the analysis selection.'
      .format(len(segsel), segsel.stop_iter-segsel.start_iter))
tts = TrajTreeSet(segsel, data_manager)
print('There are {:d} segments in the tree map.'.format(len(tts)))

tts.trace_trajectories(fppw.visit)
print('First formation counts:', fppw.n_first)
print('Unnormalized/normalized probability of first formation:', fppw.p_first, fppw.p_first/fppw.p_first.sum())
            
        

