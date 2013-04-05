# cython: profile=False
from __future__ import print_function,division
import cython
import numpy
cimport numpy

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t
ctypedef numpy.int64_t seg_id_t

weight_dtype = numpy.float64  
index_dtype = numpy.uint16
bool_dtype = numpy.bool_

from westpa.binning.assign import UNKNOWN_INDEX as _UNKNOWN_INDEX
cdef index_t UNKNOWN_INDEX = _UNKNOWN_INDEX  

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef flux_assign(numpy.ndarray[weight_t, ndim=1] weights,
                  numpy.ndarray[index_t, ndim=1] init_assignments,
                  numpy.ndarray[index_t, ndim=1] final_assignments,
                  numpy.ndarray[weight_t, ndim=2] flux_matrix):
    cdef:
        Py_ssize_t m,n
        index_t i, j
    n = len(weights)
    for m from 0 <= m < n:
        i = init_assignments[m]
        j = final_assignments[m]
        flux_matrix[i,j] += weights[m]
    return
                
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef pop_assign(numpy.ndarray[weight_t, ndim=1] weights,
                 numpy.ndarray[index_t, ndim=1] assignments,
                 numpy.ndarray[weight_t, ndim=1] populations):
    cdef:
        Py_ssize_t m,n
        index_t i,
    n = len(weights)
    for m from 0 <= m < n:
        i = assignments[m]
        populations[i] += weights[m]
    return

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef calc_rates(numpy.ndarray[weight_t, ndim=3] fluxes,
                 numpy.ndarray[weight_t, ndim=2] populations,
                 numpy.ndarray[weight_t, ndim=3] rates,
                 numpy.ndarray[bool_t,ndim=2,cast=True] masks):
    '''Calculate a series of rate matrices from a series of flux and population matrices.
    A set of boolean matrices, of the same shape as populations, is also produced, to be 
    used for generating masks for the rate matrices where initial state populations are
    zero.'''
    
    cdef:
        Py_ssize_t narrays, nbins
        index_t iarray, i, j
        
    narrays = fluxes.shape[0]
    nbins = fluxes.shape[1]
    
    for iarray from 0 <= iarray < narrays:
        for i from 0 <= i < nbins:
            if populations[iarray,i] == 0.0:
                masks[iarray,i] = 0
                for j from 0 <= j < nbins:
                    rates[iarray,i,j] = 0
            else:
                masks[iarray,i] = 1
                for j from 0 <= j < nbins:
                    rates[iarray,i,j] = fluxes[iarray,i,j] / populations[iarray,i]
    return

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef weight_t calculate_labeled_fluxes_alllags(Py_ssize_t nstates, 
                                                weights,
                                                parent_ids,
                                                micro_assignments,
                                                traj_assignments,
                                                weight_t[:,:,:,:] fluxes) except 0.0:
    cdef:
        Py_ssize_t niters = len(weights), nsegs, npts
        weight_t twindow = 0.0
        long lastiter, firstiter, iiter, windowlen
        long seg_id, current_id, parent_id
        index_t ibin, ilabel, fbin, flabel
        weight_t weight
        weight_t[:] lweights
        
        index_t[:,:] lmicro
        index_t[:,:] ltraj
    
    # We need to trace backward in each window, so we go from end to beginning
    
    for lastiter in range(niters-1,-1,-1):
        for windowlen in range(1,niters+1):
            firstiter = lastiter-windowlen+1
            if firstiter < 0: continue
            
            # we loop over all trajectories that are alive as of the last iteration
            # in the averaging window
            lweights = weights[lastiter]
            lmicro = micro_assignments[lastiter]
            ltraj  = traj_assignments[lastiter]
            nsegs = lmicro.shape[0]
            npts =  lmicro.shape[1]
            
            for seg_id in range(nsegs):
                weight = lweights[seg_id]
                fbin = lmicro[seg_id,npts-1]
                flabel = ltraj[seg_id,npts-1]
                
                # trace upwards in history to firstiter
                iiter = lastiter
                current_id = seg_id
                parent_id = parent_ids[iiter][seg_id]
                while iiter > firstiter and parent_id >= 0:
                    iiter -= 1
                    current_id = parent_id
                    parent_id = parent_ids[iiter][current_id]
                                
                ibin = micro_assignments[iiter][current_id][0]
                ilabel = traj_assignments[iiter][current_id][0]
                
                if ilabel >= nstates or flabel >= nstates:
                    raise ValueError('invalid state index (ilabel={},flabel={})'.format(ilabel,flabel))

                fluxes[ilabel,flabel,ibin,fbin] += weight
                twindow += weight*windowlen
    return twindow

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef weight_t calculate_labeled_fluxes(Py_ssize_t nstates, 
                                        weights,
                                        parent_ids,
                                        micro_assignments,
                                        traj_assignments,
                                        weight_t[:,:,:,:] fluxes) except 0.0:
    cdef:
        Py_ssize_t niters = len(weights), nsegs, npts
        weight_t twindow = 0.0
        long lastiter, firstiter, iiter, windowlen
        long seg_id, current_id, parent_id
        index_t ibin, ilabel, fbin, flabel
        weight_t weight
        weight_t[:] lweights
        
        index_t[:,:] lmicro
        index_t[:,:] ltraj
    
    # we loop over all trajectories that are alive as of the last iteration
    # in the averaging window
    
    lastiter = niters-1
    windowlen = niters
    
    lweights = weights[lastiter]
    lmicro = micro_assignments[lastiter]
    ltraj  = traj_assignments[lastiter]
    nsegs = lmicro.shape[0]
    npts =  lmicro.shape[1]
    
    for seg_id in range(nsegs):
        weight = lweights[seg_id]
        fbin = lmicro[seg_id,npts-1]
        flabel = ltraj[seg_id,npts-1]
        
        # trace upwards in history to firstiter
        iiter = lastiter
        current_id = seg_id
        parent_id = parent_ids[iiter][seg_id]
        while iiter > 0 and parent_id >= 0:
            iiter -= 1
            current_id = parent_id
            parent_id = parent_ids[iiter][current_id]
        
        assert iiter == 0 or parent_id < 0
        assert 0 <= iiter < niters
        assert current_id >= 0
        
        ibin = micro_assignments[iiter][current_id][0]
        ilabel = traj_assignments[iiter][current_id][0]

        if ilabel >= nstates or flabel >= nstates:
            raise ValueError('invalid state index (ilabel={},flabel={})'.format(ilabel,flabel))
                        
        fluxes[ilabel,flabel,ibin,fbin] += weight
        twindow += weight*windowlen
    return twindow


@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef object nested_to_flat_matrix(weight_t[:,:,:,:] input):
    '''Convert nested flux/rate matrix into a flat supermatrix.'''
    
    cdef:
        Py_ssize_t nstates = input.shape[0], nbins=input.shape[3], istate, ibin, jstate, jbin
        weight_t[:,:] _output
        
    output = numpy.empty((nstates*nbins,nstates*nbins), weight_dtype)
    _output = output
    
    for istate in range(nstates):
        for jstate in range(nstates):
            for ibin in range(nbins):
                for jbin in range(nbins):
                    #_output[istate*nbins+ibin, jstate*nbins+jbin] = input[istate, jstate, ibin, jbin]
                    _output[ibin*nstates+istate,jbin*nstates+jstate] = input[istate,jstate,ibin,jbin]
    
    return output

@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef object nested_to_flat_vector(weight_t[:,:] input):
    '''Convert nested labeled population vector into a flat vector.'''
    
    cdef:
        Py_ssize_t nstates = input.shape[0], nbins=input.shape[1], istate, ibin
        weight_t[:] _output
        
    output = numpy.empty((nstates*nbins,), weight_dtype)
    _output = output
    
    for istate in range(nstates):
        for ibin in range(nbins):
            #_output[istate*nbins+ibin] = input[istate, ibin]
            _output[ibin*nstates+istate] = input[istate,ibin]
            
    return output

@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef object flat_to_nested_matrix(Py_ssize_t nstates, Py_ssize_t nbins, weight_t[:,:] input):
    '''Convert flat supermatrix into nested matrix.'''
    
    cdef:
        Py_ssize_t istate, jstate, ibin, jbin
        weight_t[:,:,:,:] _output
        
    if input.shape[0] != nstates*nbins or input.shape[1] != nstates*nbins:
        # since input.shape is a C vector rather than a tuple, we can't print
        # it easily
        raise TypeError('input has incorrect shape for {} states and {} bins'.format(nstates, nbins))
    
    output = numpy.empty((nstates, nstates, nbins, nbins), weight_dtype)
    _output = output
    
    for istate in range(nstates):
        for jstate in range(nstates):
            for ibin in range(nbins):
                for jbin in range(nbins):
                    #_output[istate,jstate,ibin,jbin] = input[istate*nbins+ibin, jstate*nbins+jbin]
                    _output[istate,jstate,ibin,jbin] = input[ibin*nstates+istate,jbin*nstates+jstate]
    return output

@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef object flat_to_nested_vector(Py_ssize_t nstates, Py_ssize_t nbins, weight_t[:] input):
    '''Convert flat "supervector" into nested vector.'''
    
    cdef:
        Py_ssize_t istate, ibin
        weight_t[:,:] _output
        
    if input.shape[0] != nstates*nbins:
        raise TypeError('input has incorrect shape for {} states and {} bins'.format(nstates, nbins))

    output = numpy.empty((nstates, nbins), weight_dtype)
    _output = output
    
    for istate in xrange(nstates):
        for ibin in xrange(nbins):
            #_output[istate,ibin] = input[istate*nbins+ibin]
            _output[istate,ibin] = input[ibin*nstates+istate]
    
    return output

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def labeled_flux_to_rate(weight_t[:,:,:,:] labeled_fluxes, weight_t[:,:] labeled_pops):
    '''Convert a labeled flux matrix and corresponding labeled bin populations to
    a labeled rate matrix.'''
    
    cdef:
        Py_ssize_t istate, jstate, ibin, jbin, nstates, nbins
        weight_t[:,:,:,:] _rates
    
    nstates = labeled_pops.shape[0]
    nbins = labeled_pops.shape[1]
    rates = numpy.empty_like(labeled_fluxes)
    _rates = rates
    
    with nogil:
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                for ibin in xrange(nbins):
                    for jbin in xrange(nbins):
                        if labeled_pops[istate,ibin] == 0.0:
                            if labeled_fluxes[istate,jstate,ibin,jbin] > 0.0:
                                with gil:
                                    raise ValueError('flux matrix nonzero but population zero')
                            else:
                                _rates[istate,jstate,ibin,jbin] = 0.0
                        else:
                            _rates[istate,jstate,ibin,jbin] = labeled_fluxes[istate,jstate,ibin,jbin] / labeled_pops[istate,ibin]
    return rates

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sequence_macro_flux_to_rate(weight_t[:,:,:] fluxes, weight_t[:,:] traj_ens_pops):
    '''Convert a sequence of macrostate fluxes and corresponding list of trajectory ensemble populations
    to a sequence of rate matrices'''
    
    cdef:
        Py_ssize_t iiter, istate, jstate, nstates
        weight_t[:,:,:] _rates
        
    rates = numpy.empty((fluxes.shape[0], fluxes.shape[1], fluxes.shape[2]), dtype=weight_dtype)
    _rates = rates
    
    with nogil:
        for iiter in xrange(fluxes.shape[0]):
            for istate in xrange(fluxes.shape[1]):
                for jstate in xrange(fluxes.shape[2]):
                    if traj_ens_pops[iiter,istate] > 0:
                        _rates[iiter,istate,jstate] = fluxes[iiter,istate,jstate] / traj_ens_pops[iiter,istate]
                    elif fluxes[iiter,istate,jstate] > 0:
                        with gil:
                            raise ValueError('flux matrix nonzero but population zero')
                    else:
                        _rates[iiter,istate,jstate] = 0
    return rates

"""
In the following ``state`` is a 4-tuple of the following arrays of doubles:
last_time[nsegs]
last_entries[nsegs,nstates]
last_exits[nsegs,nstates]
last_completions[nsegs,nstates]


It is intended to be opaque the the calling routines
"""

@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef _fast_transition_state_copy(Py_ssize_t iiter,
                                  Py_ssize_t nstates, 
                                  seg_id_t[:] parent_ids,
                                  object last_state):
    cdef:
        bint has_last_state = 0
        Py_ssize_t nsegs, seg_id, parent_id
        double[:] _last_time, _prev_last_time
        double[:,:] _last_entries, _last_exits, _prev_last_entries, _prev_last_exits
        double[:,:,:] _last_completions, _prev_last_completions
        
    
    nsegs = parent_ids.shape[0]
    
    last_time = numpy.empty((nsegs,), numpy.double)
    last_entries = numpy.empty((nsegs,nstates), numpy.double)
    last_exits = numpy.empty((nsegs,nstates), numpy.double)
    last_completions = numpy.empty((nsegs,nstates,nstates), numpy.double)
    
    _last_time = last_time
    _last_entries = last_entries
    _last_exits = last_exits
    _last_completions = last_completions
    
    has_last_state = (last_state is not None)
    
    if has_last_state:
        _prev_last_time = last_state[0]
        _prev_last_entries = last_state[1]
        _prev_last_exits = last_state[2]
        _prev_last_completions = last_state[3]
    
    for seg_id in xrange(nsegs):
        parent_id = parent_ids[seg_id]
                
        if not has_last_state or parent_id < 0:
            _last_time[seg_id] = 0.0
            _last_entries[seg_id,:] = 0.0
            _last_exits[seg_id,:] = 0.0
            _last_completions[seg_id,:,:] = 0.0
        else:
            _last_time[seg_id] = _prev_last_time[parent_id]
            _last_entries[seg_id,:] = _prev_last_entries[parent_id,:]
            _last_exits[seg_id,:] = _prev_last_exits[parent_id,:]
            _last_completions[seg_id,:,:] = _prev_last_completions[parent_id,:,:]
            
    return (last_time, last_entries, last_exits, last_completions)


@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef find_macrostate_transitions(Py_ssize_t nstates, 
                                  weight_t[:] weights,
                                  index_t[:,:] label_assignments, 
                                  double dt, 
                                  object state,
                                  weight_t[:,:] macro_fluxes, 
                                  object durations):
    cdef:
        Py_ssize_t nsegs, npts, seg_id, ipt
        double itime, tm, t_ed
        double[:] _last_time
        double[:,:] _last_entries, _last_exits
        double[:,:,:] _last_completions
        index_t flabel, ilabel, iistate
        weight_t _weight
        

    nsegs = label_assignments.shape[0]
    npts = label_assignments.shape[1]
    
    _last_time = state[0]
    _last_entries = state[1]
    _last_exits = state[2]
    _last_completions = state[3]
    
    for seg_id in xrange(nsegs):
        itime = _last_time[seg_id]
        _weight = weights[seg_id]

        # transitions never occur between the (overlapping) end point of previous iteration and beginning of
        # current iteration, so it suffices to start looking at timepoint 1 (and backwards to timepoint 0)
        for ipt in range(1,npts):
            tm = itime + ipt*dt
            flabel = label_assignments[seg_id,ipt]
            ilabel = label_assignments[seg_id,ipt-1]
            
            # if we have wound up in a new kinetic macrostate
            if flabel != ilabel:
                for iistate in xrange(nstates):
                    if iistate != flabel:
                        macro_fluxes[iistate,flabel] += _weight
                        t_ed = tm - _last_exits[seg_id,iistate]
                        durations.append((iistate,flabel,t_ed,_weight))
                    _last_completions[seg_id,iistate,flabel] = tm
                _last_exits[seg_id,ilabel] = tm
                _last_entries[seg_id,flabel] = tm
        
        _last_time[seg_id] = tm
        

