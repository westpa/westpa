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

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef accumulate_labeled_populations(weight_t[:]  weights,
                                     index_t[:,:] bin_assignments,
                                     index_t[:,:] label_assignments,
                                     weight_t[:]  bin_pops,
                                     weight_t[:,:] labeled_bin_pops):
    '''For a set of segments in one iteration, calculate the average population in each bin, with and without
    separation by last-visited macrostate.'''
    cdef:
        Py_ssize_t nsegs, npts, nstates, seg_id, ipt
        index_t assignment, traj_assignment
        weight_t ptwt
    
    nsegs = bin_assignments.shape[0]
    npts = bin_assignments.shape[1]
    nstates = labeled_bin_pops.shape[0]
    
    with nogil:
        for seg_id in range(nsegs):
            ptwt = weights[seg_id] / npts
            for ipt in range(npts):
                assignment = bin_assignments[seg_id,ipt]
                if assignment == UNKNOWN_INDEX:
                    with gil:
                        raise ValueError('invalid bin assignment for segment {} point {}'.format(seg_id, ipt))
                else:
                    bin_pops[assignment] += ptwt
                    traj_assignment = label_assignments[seg_id,ipt]
                    if traj_assignment != UNKNOWN_INDEX:
                        labeled_bin_pops[traj_assignment,assignment] += ptwt

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
                                                weight_t[:,:,:,:] fluxes):
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
                
                assert iiter == 0 or parent_id < 0
                assert 0 <= iiter < niters
                assert current_id >= 0
                
                ibin = micro_assignments[iiter][current_id][0]
                ilabel = traj_assignments[iiter][current_id][0]
                                
                #fluxes[ilabel*nstates+ibin,flabel*nstates+fbin] += weight
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
                                        weight_t[:,:,:,:] fluxes):
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
                        
        #fluxes[ilabel*nstates+ibin,flabel*nstates+fbin] += weight
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
