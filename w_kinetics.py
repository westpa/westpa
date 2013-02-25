'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from copy import copy, deepcopy
from west.data_manager import seg_id_dtype
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection, BinMappingComponent
from collections import deque
import numpy, h5py
import scipy.linalg
import scipy.sparse.linalg
from h5py import h5s
from westtools import h5io

from west.binning.assign import index_dtype



import numba
from numba import f8, i8 #@UnresolvedImport
from pandas.io.data import DataReader

log = logging.getLogger('westtools.w_kinetics')

weight_dtype = numpy.float64
timepoint_dtype = numpy.int64

def flux_to_tprob(fluxes, pops):
    '''Calculate the transition probability matrix from flux and population'''
    tprob = numpy.empty_like(fluxes)
    for i in xrange(len(fluxes)):
        tprob[i] = fluxes[i] / pops[i]
        tprob[i] /= tprob[i].sum()
    return tprob

def get_steady_state(T):
    T = numpy.asmatrix(T)
    #ivec = numpy.empty((T.shape[0],), numpy.float64)
    #ivec[:] = 1.0/T.shape[0]
    #return (T**10000)*numpy.matrix(ivec).T
    ss = numpy.diag(T**10000).copy()
    ss /= ss.sum()
    return ss

def get_steady_state_eig(T):
    try:
        _, vecs = scipy.sparse.linalg.eigs(T.T,k=1,which='LM')
    except scipy.sparse.linalg.eigen.arpack.ArpackError:
        return get_steady_state(T)
    else:
        ss = numpy.abs(vecs[:,0])
        ss /= ss.sum()
        return ss

@numba.jit(restype=numba.object_,
           argtypes=[f8,f8,f8,numba.uint16,
                     numba.uint8[:], numba.uint16[:],
                     f8[:], f8[:],
                     f8[:,:],f8[:,:],
                     f8[:],f8[:],f8[:,:]],
           locals=dict(new_last_exit=f8[:], new_last_entry=f8[:], new_last_completion=f8[:,:]))

def find_transitions_mz(weight, initial_time, dt, last_macro,
                        macromask, bins, 
                        pops, macro_pops, 
                        micro_fluxes, macro_fluxes,
                        last_exit, last_entry, last_completion):
    '''Find transitions, updating population and flux data in place and returning a reference to the new 
    per-trajectory state dict.'''
    
    npts  = len(bins)
    nbins = len(pops)
    popwt = 1.0/(npts-1)
    
    # Preserve old history
    new_last_exit = last_exit.copy()
    new_last_entry = last_entry.copy()
    new_last_completion = last_completion.copy()

    event_durations = []
    
    # transitions never occur during the (overlapping) end point of previous iteration and
    # beginning point of current iteration, so start scan at point 1 instead of point 0
    for ipt in xrange(1,npts):
        tm = initial_time + dt*ipt # current time
        
        fbin = bins[ipt]
        ibin = bins[ipt-1]
        
        pops[ibin] += weight*popwt
        
        # We always update the microstate flux matrix
        micro_fluxes[ibin,fbin] += weight
                
        if fbin != ibin:
            # transition has occurred; record times, etc
            new_last_exit[ibin] = tm
            new_last_entry[fbin] = tm
            
            # we only count transitions ending in a kinetic macrostate
            if macromask[fbin] == 1:
                # loop over all possible initial states
                for iibin in xrange(nbins):
                    if (# we only count transitions from non-transition regions
                        macromask[iibin] == 1 
                        
                        # and only those where walkers actually departed from at some point
                        and new_last_exit[iibin] > 0 
                        
                        # and we need to have visited the initial state more recently than this final state,
                        # or else we would be double-counting transitions
                        and new_last_entry[iibin] > new_last_completion[iibin,fbin] 
                        ):
                        # exclude sojourns out of bin and back for now, because it clutters up
                        # output. TODO: this should be an option
                        if iibin != fbin:
                            macro_fluxes[iibin,fbin] += weight
                            t_ed = new_last_entry[fbin] - new_last_exit[iibin]
                            event_durations.append((iibin,fbin,t_ed,weight))
                        new_last_completion[iibin,fbin] = tm
                    last_macro = fbin
                    
    if macromask[last_macro]:
        macro_pops[last_macro] += weight
        
    return (event_durations,
            {'last_macrostate': last_macro,
             'last_entry': new_last_entry,
             'last_exit': new_last_exit, 
             'last_completion': new_last_completion,
             'initial_time': tm}
            )

@numba.jit(argtypes=[numba.uint16[:,:], f8[:], f8[:]])
def pops_from_assignments(assignments, weights, pops):
    nsegs, npts = assignments.shape
    
    for seg_id in xrange(nsegs):
        tweight = weights[seg_id]/npts 
        for ipt in xrange(npts):
            pops[assignments[seg_id,ipt]] += tweight
    
#@profile
def build_flux_matrix(data_reader, assignments_file, max_history):
    max_history = 5
    
    nbins = assignments_file.attrs['nbins']
    start_iter, stop_iter = h5io.get_iter_range(assignments_file)
    iter_count = stop_iter - start_iter
    
    pops_by_iter = numpy.zeros((iter_count,nbins), weight_dtype)
    rate_matrix_by_iter = numpy.zeros((iter_count,nbins,nbins), weight_dtype)
    flux_matrix_by_iter = numpy.zeros((iter_count,nbins,nbins), weight_dtype)
    
    last_history = None
    for iiter,n_iter in enumerate(xrange(start_iter,100)):#stop_iter)):
        print(iiter,n_iter)
        
        iter_group = data_reader.get_iter_group(n_iter)
        seg_index = iter_group['seg_index'] # weight, parent_id
        weights = seg_index['weight']
        parent_ids = seg_index['parent_id']            
        n_segs = weights.shape[0]
        
        assignments_ds = assignments_file.get_iter_group(n_iter)['assignments']
        assignments_data = assignments_ds[...]
        
        pops_from_assignments(assignments_data, weights, pops_by_iter[iiter])
        
        next_history = [None]*n_segs
        windowlen = 0.0 # number of tau accumulated, can be fractional if some trajectories don't cover entire window
        for seg_id in xrange(n_segs):
            parent_id = parent_ids[seg_id]
            weight = weights[seg_id]
            
            if parent_ids[seg_id] < 0:
                traj_history = deque(maxlen=max_history)
            else:
                traj_history = copy(last_history[parent_id])
                
            traj_history.appendleft(assignments_data[seg_id,[0,-1]])            
                                    
            windowlen += _flux_acc(traj_history, max_history, flux_matrix_by_iter[iiter], weight)
            next_history[seg_id] = traj_history
        
        n = max_history
        timefac = (n+1)*(n+1)*n/2 - n*(n+1)*(2*n+1)/6
        #print(windowlen, timefac)
        flux_matrix_by_iter[iiter] /= windowlen
        print(pops_by_iter[iiter])
        print(flux_matrix_by_iter[iiter])
        
        rate_matrix_by_iter[iiter] = flux_matrix_by_iter[iiter]
        
        for i in xrange(nbins):
            if pops_by_iter[iiter,i] == 0: continue
            print(pops_by_iter[iiter,i])
            rate_matrix_by_iter[iiter,i] /= pops_by_iter[iiter,i]
            
        print(rate_matrix_by_iter[iiter])
            
        
        # prepare for next iteration
        last_history = next_history
        del assignments_data
    
    avg_rate_matrix = rate_matrix_by_iter[:98].mean(axis=0)
    T = avg_rate_matrix.copy()
    avg_pops = pops_by_iter[:98].mean(axis=0)
    
    print('avg rate matrix')
    print(avg_rate_matrix)
    
    for i in xrange(nbins):
        T[i] /= T[i].sum()
        
    ss_pops = get_steady_state_eig(T)
    print('iteration', n_iter)
    print('avg pops')
    print(avg_pops)
    print('ss pops')
    print(ss_pops)
    
    # A = unbound = 2
    # B = bound = 0 
    
    kAB = 0.0
    kBA = 0.0
    for i in xrange(nbins):
        kAB += avg_rate_matrix[i,0] * ss_pops[i]
        kBA += avg_rate_matrix[i,2] * ss_pops[i]
        
    print('kAB', kAB)
    print('kBA', kBA)
    
    

"""
@numba.jit(restype=f8,
           argtypes=[numba.object_, i8, f8[:,:], f8],
           locals=dict(acc_tau=f8, hist_len=i8))
"""
def _flux_acc(traj_history, max_history, flux_matrix, weight):
    acc_tau = 0
    hist_len = int(len(traj_history))
    
    # so much potential for off-by-one errors here
    # don't mess with this without testing            
        
    # for windows of length [1,minimum of available iterations, max history length]            
    for windowlen in xrange(1,hist_len+1):
        # for offset within history window [0,max history length-1]
        for offset in xrange(max_history):
            if offset+windowlen-1 > hist_len-1: continue
                                    
            ibin = traj_history[windowlen+offset-1][1]
            fbin = traj_history[offset][0]
            flux_matrix[ibin,fbin] += weight
            
            acc_tau += windowlen*weight
    return acc_tau


class WKinetics(WESTTool):
    prog='w_kinetics'
    description = '''\
Calculate kinetics data for arbitrary states from weighted ensemble data.

A bin assignment file (usually "assign.h5") is required. This file should
be produced with w_assign such that each kinetic macrostate corresponds
to exactly one bin, and those bins should have been identified as kinetic
macrostates (see "w_assign --help" for information).

'''
    
    def __init__(self):
        super(WKinetics,self).__init__()
        self.data_reader = WESTDataReader() 
        self.output_file = None
        self.assignments_file = None
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        parser.add_argument('-a', '--assignments', default='assign.h5',
                            help='''Bin assignments and macrostate definitions are in ASSIGNMENTS
                            (default: %(default)s).''')
        parser.add_argument('-o', '--output', dest='output', default='kinetics.h5',
                            help='''Store results in OUTPUT (default: %(default)s).''')

        
    def process_args(self, args):
        self.assignments_file = h5io.WESTPAH5File(args.assignments, 'r')
        self.data_reader.process_args(args)
        self.data_reader.open('r')
        self.output_file = h5io.WESTPAH5File(args.output, 'w', creating_program=True)
        h5io.stamp_creator_data(self.output_file)

    def walk_tree_mz(self):        

        current_state = None # list of dicts
        nbins = self.assignments_file.attrs['nbins']
        start_iter, stop_iter = h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
        
        avg_micro_flux = numpy.zeros((nbins,nbins), dtype=weight_dtype)
        avg_macro_flux = numpy.zeros((nbins,nbins), dtype=weight_dtype)
        avg_pops = numpy.zeros((nbins,), dtype=weight_dtype)
        all_eds = []
        
        avg_macro_pops = numpy.zeros((nbins,), dtype=weight_dtype)
        
        macromask = numpy.zeros((nbins,), dtype=numpy.bool_)
        for macrostate_bin in self.assignments_file['state_assignments']:
            macromask[macrostate_bin] = 1
        
        for iiter,n_iter in enumerate(xrange(start_iter,stop_iter)):
            print(iiter,n_iter)
            
            iter_group = self.data_reader.get_iter_group(n_iter)
            seg_index = iter_group['seg_index'] # weight, parent_id
            weights = seg_index['weight']
            parent_ids = seg_index['parent_id']            
            pcoord_ds = iter_group['pcoord']
            n_segs, n_points = pcoord_ds.shape[:2]
            assignment_data = self.assignments_file.get_iter_group(n_iter)['assignments'][...]
            
            dt = 1.0/(n_points-1)
            
            pops = numpy.zeros((nbins,), dtype=weight_dtype)
            macro_pops = numpy.zeros((nbins,), dtype=weight_dtype)
            micro_fluxes = numpy.zeros((nbins,nbins), dtype=weight_dtype)
            macro_fluxes = numpy.zeros((nbins,nbins), dtype=weight_dtype)
            
            next_state = [None] * n_segs
            for seg_id in xrange(n_segs):
                weight = weights[seg_id]
                parent_id = parent_ids[seg_id]
                bins = assignment_data[seg_id]
                                    
                if parent_id < 0:
                    state = {'last_macrostate': bins[0],
                             'last_exit': numpy.zeros((nbins,), dtype=numpy.float64),
                             'last_entry': numpy.zeros((nbins,), dtype=numpy.float64),
                             'last_completion': numpy.zeros((nbins,nbins), dtype=numpy.float64),
                             'initial_time': 0.0}
                else:
                    assert current_state[parent_id] is not None
                    state = current_state[parent_id]
                    
                eds, state = find_transitions_mz(weight, state['initial_time'], dt, state['last_macrostate'],
                                                 macromask, bins, 
                                                 pops, macro_pops,
                                                 micro_fluxes, macro_fluxes, 
                                                 state['last_exit'], state['last_entry'], state['last_completion'])
                all_eds.extend(eds)
                if eds:
                    print(eds)
                next_state[seg_id] = state
                
                del state, bins
            current_state = next_state
            
            avg_micro_flux += micro_fluxes
            avg_macro_flux += macro_fluxes
            avg_pops += pops 
            avg_macro_pops += macro_pops/macro_pops.sum()
            
            # clean up for next iteration
            del seg_index, weights, parent_ids, pcoord_ds, micro_fluxes, macro_fluxes, assignment_data

        print('number of durations recorded:',len(all_eds))

        avg_pops /= iter_count
        avg_micro_flux /= iter_count
        avg_macro_pops /= iter_count
        avg_macro_flux /= iter_count

        print('\npops:')
        print(repr(avg_pops))
        print('norm check:', avg_pops.sum())
        
        print('\nmacropops:')
        print(repr(avg_macro_pops))
        print('norm check:', avg_macro_pops.sum())
        
        print('\nmicrostate (bin-to-bin) flux matrix:')
        print(repr(avg_micro_flux))
        
        tprob = flux_to_tprob(avg_micro_flux, avg_pops)
        if not numpy.isnan(tprob).any():
            print('microstate (bin-to-bin) transition probability matrix:')
            print(repr(tprob))
            try:
                #vals, vects = scipy.linalg.eig(tprob,left=True,right=False)
                vals, vects = scipy.sparse.linalg.eigs(tprob.T,k=1,which='LM')
            except ValueError:
                print('no steady state solution')
            else:                
                ss = get_steady_state(tprob)
                print('steady-state solution (from powering):',ss)
                print('norm check:', ss.sum())
                asorted = numpy.argsort(numpy.abs(vals))
                ss = numpy.abs(vects[:,asorted[-1]])
                print('steady-state solution (from eigenvalues):',ss)
                print('norm check:', ss.sum())
                print('renormalized from eigenvalues:',ss/ss.sum())
                
                ss/=ss.sum()
                print('norm of difference between avg and ss pop vectors:', (((pops-ss)**2).sum())**0.5)
                
            
                
        print('\nmacrostate flux matrix:')
        print(repr(avg_macro_flux))
        
    def build_rate_matrix(self):
        build_flux_matrix(self.data_reader, self.assignments_file, 50)
                    
            
            
            
            
        
        
    def go(self):
        self.build_rate_matrix()
        #self.walk_tree_mz()
        

if __name__ == '__main__':
    WKinetics().main()
