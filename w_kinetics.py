'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection
from collections import namedtuple, deque
import numpy, h5py
import scipy.linalg
import scipy.sparse.linalg
from westtools import h5io

from west.binning.assign import index_dtype, UNKNOWN_INDEX
from west.data_manager import seg_id_dtype, n_iter_dtype, weight_dtype

import numba
from numba import f8, i8, u8 #@UnresolvedImport

log = logging.getLogger('westtools.w_kinetics')

_rate_result = namedtuple('_rate_result', ('micro_pops','traj_pops','micro_fluxes', 'micro_rates', 'macro_rates'))

@numba.jit(argtypes=[f8[:], numba.uint16[:,:], numba.uint16[:,:], f8[:], f8[:]])
def _acc_pops(weights, micro_assignments, traj_assignments, micro_pops, traj_pops):
    nsegs = micro_assignments.shape[0]
    npts  = micro_assignments.shape[1]
    for seg_id in xrange(nsegs):
        for ipt in xrange(npts):
            micro_pops[micro_assignments[seg_id,ipt]] += weights[seg_id] / npts
            # again, account for possibly broken numba where 65535 -> -1 when it shouldn't
            if traj_assignments[seg_id,ipt] != UNKNOWN_INDEX and traj_assignments[seg_id,ipt] != -1:
                traj_pops[traj_assignments[seg_id,ipt]] += weights[seg_id] / npts

@numba.jit(argtypes=[i8, numba.object_, numba.object_, numba.object_, numba.object_, numba.object_],
           locals=dict(weight=f8, ibin=numba.uint16, ilabel=numba.uint16, fbin=numba.uint16, flabel=numba.uint16,
                       iiter=i8, firstiter=i8, lastiter=i8, parent_id=i8, current_id=i8, windowlen=i8, twindow=f8, fluxes=f8[:,:]))                
def estimate_rates(nbins, state_labels, weights, parent_ids, micro_assignments, traj_assignments):
    '''Estimate rates over multiple iterations. All possible windows are used.
    Returns microstate (bin) populations, trajectory (last-in-state) populations, 
    microstate (bin-to-bin) flux matrix, microstate rate matrix,
    macrostate (state-to-state) rate matrix.'''
    
    assert len(weights) == len(parent_ids) == len(micro_assignments) == len(traj_assignments)
    
    niters = len(weights)
    nstates = len(state_labels)
        
    # Estimate microstate and trajectory populations over entire window
    micro_pops = numpy.zeros((nbins,), weight_dtype)
    traj_pops =  numpy.zeros((nstates,), weight_dtype)
    
    for iiter in xrange(niters):
        _acc_pops(weights[iiter], micro_assignments[iiter], traj_assignments[iiter], micro_pops, traj_pops)
        
    micro_pops /= niters
    traj_pops /= niters 
    
    # Loop over all possible windows to accumulate flux matrix
    fluxes = numpy.zeros((nbins*nstates,nbins*nstates), weight_dtype)
    
    # We need to trace backward in each window, so we go from end to beginning
    twindow = 0.0
    for lastiter in xrange(niters-1,-1,-1):
        for windowlen in xrange(1,niters+1):
            firstiter = lastiter-windowlen+1
            if firstiter < 0: continue
            
            # we loop over all trajectories that are alive as of the last iteration
            # in the averaging window
            nsegs = len(weights[lastiter])
            for seg_id in xrange(nsegs):
                weight = weights[lastiter][seg_id]
                fbin = micro_assignments[lastiter][seg_id][-1]
                flabel = traj_assignments[lastiter][seg_id][-1]
                
                # trace upwards in history to firstiter
                iiter = lastiter
                current_id = seg_id
                parent_id = parent_ids[iiter][seg_id]
                while iiter > firstiter and parent_id >= 0:
                    iiter -= 1
                    current_id = parent_id
                    parent_id = parent_ids[iiter][current_id]
                assert iiter == firstiter or parent_id < 0
                
                ibin = micro_assignments[iiter][current_id][0]
                ilabel = traj_assignments[iiter][current_id][0]
                                
                fluxes[nstates*ibin+ilabel,nstates*fbin+flabel] += weight
                twindow += weight*windowlen
    log.debug('averaging window = {} tau'.format(twindow))

    
    fluxes /= twindow    
    print(fluxes)

class WKinetics(WESTTool):
    prog='w_kinetics'
    description = '''\
Calculate kinetics data for arbitrary states from weighted ensemble data.

A bin assignment file (usually "assign.h5") including trajectory labeling
is required (see "w_assign --help" for information on generating this file).

Flux and rate matrices, along with populations, may be estimated over
windows of several tau. This averaging is performed in such a way that
all possible windows of that size and smaller are used. (For instance, with a
window of 5 iterations, the rate matrix is built using 5 1-iteration windows,
4 2-iteration windows, 3 3-iteration windows, 2 4-iteration windows, and 1 
5-iteration window ending at each tau.)
'''
    
    def __init__(self):
        super(WKinetics,self).__init__()
        self.data_reader = WESTDataReader() 
        self.output_file = None
        self.assignments_file = None
        self.window_size = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        parser.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')
        parser.add_argument('-o', '--output', dest='output', default='kinetics.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')
        parser.add_argument('-w', '--windowsize', type=int, default=1,
                            help='''Estimate kinetics over a maximum of WINDOWSIZE iterations.
                            (Default: %(default)s).''')

        
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.data_reader.process_args(args)
        self.data_reader.open('r')
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)
        self.window_size = args.windowsize
        
    def go(self):
        nbins = self.assignments_file.attrs['nbins']
        state_labels = self.assignments_file['state_labels']
        start_iter, stop_iter = h5io.get_iter_range(self.assignments_file)
        
        weights_ring = deque(maxlen=self.window_size)
        parent_ids_ring = deque(maxlen=self.window_size)
        micro_assignments_ring = deque(maxlen=self.window_size)
        traj_assignments_ring = deque(maxlen=self.window_size)
        
        for iiter, n_iter in enumerate(xrange(start_iter,10)):# stop_iter)):
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index']
            nsegs, npts = iter_group['pcoord'].shape[0:2] 
            weights = seg_index['weight']
            parent_ids = seg_index['parent_id']
            micro_assignments = self.assignments_file['assignments'][iiter,:nsegs,:npts]
            traj_assignments = self.assignments_file['trajlabels'][iiter,:nsegs,:npts]
            
            weights_ring.append(weights)
            parent_ids_ring.append(parent_ids)
            micro_assignments_ring.append(micro_assignments)
            traj_assignments_ring.append(traj_assignments)
            
            print(n_iter)
            estimate_rates(nbins, state_labels, weights_ring, parent_ids_ring, micro_assignments_ring, traj_assignments_ring)
            
            del iter_group, weights, parent_ids, micro_assignments, traj_assignments
            
        

if __name__ == '__main__':
    WKinetics().main()
