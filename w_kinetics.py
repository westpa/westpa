'''
Created on Feb 15, 2013

@author: mzwier
'''

from __future__ import print_function, division; __metaclass__ = type
import logging
from copy import copy, deepcopy
from west.data_manager import seg_id_dtype
from westtools.tool_classes import WESTTool, WESTDataReader, IterRangeSelection, BinMappingComponent
import numpy, h5py
import scipy.linalg
from h5py import h5s
from westtools import h5io

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
    return numpy.diag(T**10000)

@numba.jit(restype=numba.object_,
           argtypes=[f8,f8,f8,numba.uint16,
                     numba.uint8[:], numba.uint16[:],
                     f8[:], f8[:],
                     f8[:,:],f8[:,:],
                     f8[:],f8[:],f8[:,:]],
           locals=dict(new_last_exit=f8[:], new_last_entry=f8[:], new_last_completion=f8[:,:]))

def find_transitions(weight, initial_time, dt, last_macro,
                     macromask, bins, 
                     pops, macro_pops, 
                     micro_fluxes, macro_fluxes,
                     last_exit, last_entry, last_completion):
    '''Find transitions, updating population and flux data in place and returning a reference to the new 
    per-trajectory state dict.'''
    
    npts  = len(bins)
    nbins = len(pops)
    popwt = 1.0/(npts-1)
        
    # We never update the history in-place, so shadow these with copies
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
            
            # are we ending in a macrostate bin?
            if macromask[fbin] == 1:
                # if so, we need to update the total flow data and calculate event duration
                # event duration is time of last exit to time of entry
                for iibin in xrange(nbins):
                    if (macromask[iibin] == 1
                        and new_last_exit[iibin] > 0 
                        and new_last_entry[iibin] > new_last_completion[iibin,fbin]
                        ):
                        macro_fluxes[iibin,fbin] += weight
                        t_ed = new_last_entry[fbin] - new_last_exit[iibin]
                        if iibin != fbin: 
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

    def walk_tree(self):        

        current_state = None # list of dicts
        nbins = self.assignments_file.attrs['nbins']
        start_iter, stop_iter = h5io.get_iter_range(self.assignments_file)
        iter_count = stop_iter - start_iter
        
        avg_micro_flux = numpy.zeros((nbins,nbins), dtype=weight_dtype)
        avg_pops = numpy.zeros((nbins,), dtype=weight_dtype)
        all_eds = []
        
        avg_macro_pops = numpy.zeros((nbins,), dtype=weight_dtype)
        
        macromask = numpy.zeros((nbins,), dtype=numpy.uint8)
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
            pcoord_data = pcoord_ds[...]
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
                    
                eds, state = find_transitions(weight, state['initial_time'], dt, state['last_macrostate'],
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
            avg_pops += pops 
            avg_macro_pops += macro_pops/macro_pops.sum()
            
            # clean up for next iteration
            del seg_index, weights, parent_ids, pcoord_ds, pcoord_data, micro_fluxes, macro_fluxes, assignment_data

        print('number of durations recorded:',len(all_eds))

        avg_pops /= iter_count
        avg_micro_flux /= iter_count
        avg_macro_pops /= iter_count

        
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
                vals, vects = scipy.linalg.eig(tprob,left=True,right=False)
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
        
    def go(self):
        self.walk_tree()
        

if __name__ == '__main__':
    WKinetics().main()
