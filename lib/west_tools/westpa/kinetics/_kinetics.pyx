# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function,division
import cython
import numpy
import warnings
cimport numpy

ctypedef numpy.uint16_t index_t
ctypedef numpy.float64_t weight_t
ctypedef numpy.uint8_t bool_t
ctypedef numpy.int64_t seg_id_t
ctypedef numpy.uint_t uint_t # 32 bits on 32-bit systems, 64 bits on 64-bit systems

cdef double NAN = numpy.nan 

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
cpdef calc_rates(weight_t[:,::1] fluxes,
                 weight_t[::1] populations,
                 weight_t[:,::1] rates,
                 bool_t[:,::1] mask):
    '''Calculate a rate matrices from flux and population matrices. A matrix of the same 
    shape as fluxes, is also produced, to be used for generating a mask for the rate 
    matrices where initial state populations are zero.'''
    
    cdef:
        Py_ssize_t narrays, nbins
        index_t iarray, i, j
        
    nbins = fluxes.shape[0]
    
    with nogil:
        for i in range(nbins):
            if populations[i] == 0.0:
                for j in range(nbins):
                    mask[i,j] = 1
                    rates[i,j] = 0.0
            else:
                for j in range(nbins):
                    mask[i,j] = 0
                    rates[i,j] = fluxes[i,j] / populations[i]


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

        #if ilabel >= nstates or flabel >= nstates:
        #    raise ValueError('invalid state index (ilabel={},flabel={})'.format(ilabel,flabel))
        if ilabel < nstates and flabel < nstates:
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

@cython.boundscheck(True)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef _reduce_labeled_rate_matrix_to_macro(Py_ssize_t nstates, Py_ssize_t nbins, weight_t[:,:] rates, weight_t[:] pops):
    '''Reduce a labeled microstate rate matrix into a macrostate rate matrix. This is
    for internal use, where the rates/pops vectors have been blocked by state.'''
    
    cdef:
        Py_ssize_t istate, jstate, ibin, jbin
        weight_t[:,:] _macro_rates
        weight_t sspop, rate_elem, traj_ens_pop
    
    macro_rates = numpy.zeros((nstates, nstates), numpy.float64)
    _macro_rates = macro_rates
    
    for istate in xrange(nstates):
        for jstate in xrange(nstates):
            for ibin in xrange(nbins):
                for jbin in xrange(nbins):
                    sspop = pops[ibin*nstates+istate]
                    rateelem = rates[ibin*nstates+istate,jbin*nstates+jstate]
                    _macro_rates[istate,jstate] += sspop*rateelem
                    
    # Normalize by total population in each trajectory ensemble
    for istate in xrange(nstates):
        #traj_ens_pop = pops[istate].sum()
        traj_ens_pop = 0
        for ibin in xrange(nbins):
            traj_ens_pop += pops[ibin*nstates+istate]
        
        for jstate in xrange(nstates):
            _macro_rates[istate, jstate] /=  traj_ens_pop
        
    return macro_rates


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef labeled_flux_to_rate(weight_t[:,:,:,:] labeled_fluxes, weight_t[:,:] labeled_pops, object output=None):
    '''Convert a labeled flux matrix and corresponding labeled bin populations to
    a labeled rate matrix.'''

    cdef:
        Py_ssize_t istate, jstate, ibin, jbin, nstates, nbins
        weight_t[:,:,:,:] _rates

    nstates = labeled_fluxes.shape[0]
    nbins = labeled_fluxes.shape[2]

    if output is None:
        output = numpy.empty_like(labeled_fluxes)
    _rates = output

    with nogil:
        for istate in xrange(nstates):
            for jstate in xrange(nstates):
                for ibin in xrange(nbins):
                    for jbin in xrange(nbins):
                        if labeled_pops[istate,ibin] == 0.0:
                            if labeled_fluxes[istate,jstate,ibin,jbin] > 0.0:
                                with gil:
                                    #raise ValueError('flux matrix entry nonzero but population zero')
                                    warnings.warn('flux matrix entry nonzero but population zero')

                            _rates[istate,jstate,ibin,jbin] = 0.0
                        else:
                            _rates[istate,jstate,ibin,jbin] = labeled_fluxes[istate,jstate,ibin,jbin] / labeled_pops[istate,ibin]
    return output

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef sequence_macro_flux_to_rate(weight_t[:] dataset, weight_t[:,:] pops, Py_ssize_t istate, Py_ssize_t jstate, bint pairwise=True, stride=None):
    '''Convert a sequence of macrostate fluxes and corresponding list of trajectory ensemble populations
    to a sequence of rate matrices.
    
    If the optional ``pairwise`` is true (the default), then rates are normalized according to the
    relative probability of the initial state among the pair of states (initial, final); this is
    probably what you want, as these rates will then depend only on the definitions of the states
    involved (and never the remaining states). Otherwise (``pairwise'' is false), the rates are
    normalized according the probability of the initial state among *all* other states.'''
    
    cdef:
        Py_ssize_t iiter, nstates, itersum
        weight_t[:] _rates, _fluxsum, _pairsum, _psum
        
    rates = numpy.zeros((dataset.shape[0]), dtype=weight_dtype)
    #rates = :W
    fluxsum = numpy.zeros((dataset.shape[0]), dtype=weight_dtype)
    psum = numpy.zeros((dataset.shape[0]), dtype=weight_dtype)
    pairsum = numpy.zeros((dataset.shape[0]), dtype=weight_dtype)
    _fluxsum = fluxsum
    _pairsum = pairsum
    _psum = psum
    _rates = rates
    
    # We want to modify this to be the SUM of fluxes up till this point, divided by the SUM of the population till then.
    with nogil:
        for iiter in xrange(dataset.shape[0]):
            if iiter == 0:
                # We need to catch if we haven't entered the istate or jstate yet.
                # Otherwise, we introduce a NaN into the calculation that we cannot recover from.
                if pairwise and (pops[0,istate] + pops[0,jstate]) != 0.0:
                    _psum[0] = pops[0, istate] / (pops[0,istate] + pops[0,jstate])
                else:
                    _psum[0] = pops[0, istate]
                _fluxsum[0] = dataset[0]
            else:
                if pairwise and (pops[iiter,istate] + pops[iiter,jstate]) != 0.0:
                    _psum[iiter] = (pops[iiter, istate] / (pops[iiter,istate] + pops[iiter,jstate])) + _psum[iiter-1]
                else:
                    _psum[iiter] = pops[iiter,istate] + _psum[iiter-1]
                _fluxsum[iiter] = dataset[iiter] + _fluxsum[iiter-1]
            if _psum[iiter] > 0 and _fluxsum[iiter] > 0:
                _rates[iiter] = _fluxsum[iiter] / _psum[iiter]
            else:
                _rates[iiter] = 0.0

    return rates[iiter]

"""
In the following ``state`` is a 5-tuple of the following arrays of doubles:
last_time[nsegs]
last_entries[nsegs,nstates]
last_exits[nsegs,nstates]
last_exits_td[nsegs,nstates]
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
        double[:,:] _last_entries, _last_exits, _prev_last_entries, _prev_last_exits, _last_exits_td, _prev_last_exits_td
        double[:,:,:] _last_completions, _prev_last_completions
        
    
    nsegs = parent_ids.shape[0]
    
    last_time = numpy.empty((nsegs,), numpy.double)
    # Use nstates + 1 to account for possible unknown states
    last_entries = numpy.empty((nsegs,nstates+1), numpy.double)
    last_exits = numpy.empty((nsegs,nstates+1), numpy.double)
    last_exits_td = numpy.empty((nsegs,nstates+1), numpy.double)
    last_completions = numpy.empty((nsegs,nstates+1,nstates+1), numpy.double)
    
    _last_time = last_time
    _last_entries = last_entries
    _last_exits = last_exits
    _last_exits_td = last_exits_td
    _last_completions = last_completions
    
    has_last_state = (last_state is not None)
    
    if has_last_state:
        _prev_last_time = last_state[0]
        _prev_last_entries = last_state[1]
        _prev_last_exits = last_state[2]
        _prev_last_exits_td = last_state[3]
        _prev_last_completions = last_state[4]
    
    for seg_id in xrange(nsegs):
        parent_id = parent_ids[seg_id]
                
        if not has_last_state or parent_id < 0:
            _last_time[seg_id] = 0.0
            _last_entries[seg_id,:] = 0.0
            _last_exits[seg_id,:] = 0.0
            _last_exits_td[seg_id,:] = 0.0
            _last_completions[seg_id,:,:] = 0.0
        else:
            _last_time[seg_id] = _prev_last_time[parent_id]
            _last_entries[seg_id,:] = _prev_last_entries[parent_id,:]
            _last_exits[seg_id,:] = _prev_last_exits[parent_id,:]
            _last_exits_td[seg_id,:] = _prev_last_exits_td[parent_id,:]
            _last_completions[seg_id,:,:] = _prev_last_completions[parent_id,:,:]
            
    return (last_time, last_entries, last_exits, last_exits_td, last_completions)


@cython.boundscheck(False)
@cython.wraparound(False)    
cpdef find_macrostate_transitions(Py_ssize_t nstates, 
                                  weight_t[:] weights,
                                  index_t[:,:] label_assignments,
                                  index_t[:,:] state_assignments,
                                  double dt, 
                                  object state,
                                  weight_t[:,:] macro_fluxes,
                                  uint_t[:,:] macro_counts,
                                  weight_t[:] target_fluxes,
                                  uint_t[:] target_counts,
                                  object durations):
    cdef:
        Py_ssize_t nsegs, npts, seg_id, ipt
        double itime, tm, t_ed
        double[:] _last_time
        double[:,:] _last_entries, _last_exits, _last_exits_td
        double[:,:,:] _last_completions
        index_t flabel, ilabel, iistate, slabel
        weight_t _weight
    """
    A cythoned function designed to track how long macrostate transitions take.  Requires the simulation
    to have already been binned and placed into macrostates, as appropriate.  Called by functions such as
    w_kinetics to generate kinetics information in the 'per-tau' format.

    Parameters
    ----------

    weights : weight_t
        Weights, typically from the main .h5 file of the simulation.  Likely called from the data reader
        of the calling function.
    label_assignments : index_t
        Macrostate label assignments, as compatible with those outputted by w_assign.  Bins are marked as
        states (or ignored) in a previous step.  Should be in the form  of a 'tag', or 'color'; in this 
        dataset, once a walker has been marked with a macrostate, it does not lose the macrostate 
        assigment, even upon leaving the appropriately defined state bin, until it enters another state bin.
    state_assignments : index_t
        Macrostate label assignments, but without any 'color' tagging.
    dt : double
        The number of timesteps of the system.
    state : object
       The output of _fast_transition_state_copy.
    macro_fluxes : weight_t
        An array that contains state to state fluxes.
    macro_counts : uint_t
        An array that contains the observed number of state to state fluxes.
    target_fluxes : weight_t
        An array that contains the fluxes into the target state from any state.
    target_counts : uint_t
        An array that contains the observed number of fluxes into the target state from any state.
    durations : list_like (object)
        A list containing the calculated duration information, including the iistate, flabel (state to state),
        event duration, the weight of all walkers involved, and all seg_ids.

    """
        

    nsegs = label_assignments.shape[0]
    npts = label_assignments.shape[1]
    
    _last_time = state[0]
    _last_entries = state[1]
    _last_exits = state[2]
    _last_exits_td = state[3]
    _last_completions = state[4]
    
    for seg_id in xrange(nsegs):
        itime = _last_time[seg_id]
        _weight = weights[seg_id]

        # transitions never occur between the (overlapping) end point of previous iteration and beginning of
        # current iteration, so it suffices to start looking at timepoint 1 (and backwards to timepoint 0)
        for ipt in range(1,npts):
            tm = itime + ipt*dt
            flabel = label_assignments[seg_id,ipt]
            ilabel = label_assignments[seg_id,ipt-1]
            slabel = state_assignments[seg_id,ipt]

            # if we have left our state transition barrier...
            if flabel == slabel:
                _last_exits_td[seg_id,flabel] = tm

            if flabel != ilabel:
                target_fluxes[flabel] += _weight
                target_counts[flabel] += 1
                _last_exits[seg_id,ilabel] = tm
                _last_entries[seg_id,flabel] = tm

                for iistate in xrange(nstates):
                    # if we have more recently returned to iistate than arrived at flabel from iistate,
                    # we note a new completed transition from iistate to flabel
                    # equality applies only for 0, which means we're counting an arrival from the
                    # state where the trajectory started
                    if _last_exits[seg_id, iistate] > 0 and _last_entries[seg_id,iistate] >= _last_completions[seg_id,iistate,flabel]:
                        macro_fluxes[iistate,flabel] += _weight
                        macro_counts[iistate,flabel] += 1
                        _last_completions[seg_id,iistate,flabel] = tm

                        # omit circular transitions (for now) because it causes the transition
                        # list to explode
                        if iistate != flabel:
                            t_ed = tm - _last_exits_td[seg_id,iistate]
                            durations.append((iistate,flabel,t_ed,_weight, seg_id))
        _last_time[seg_id] = tm


cdef class StreamingStats2D:
    '''Calculate mean and variance of a series of two-dimensional arrays of shape (nbins, nbins)
    using an online algorithm. The statistics are accumulated along what would be axis=0 if the 
    input arrays were stacked vertically. 

    This code has been adapted from:
    http://www.johndcook.com/skewness_kurtosis.html'''

    cdef weight_t[:,::1] _M1
    cdef weight_t[:,::1] _M2
    cdef uint_t[:,::1] _n
    cdef Py_ssize_t _sz0, _sz1

    def __init__(self, tuple shape):

        assert len(shape) == 2

        self._n = numpy.zeros(shape, dtype=numpy.uint)
        self._M1 = numpy.zeros(shape, dtype=weight_dtype)
        self._M2 = numpy.zeros(shape, dtype=weight_dtype)
        self._sz0, self._sz1 = shape

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def update(self, weight_t[:,::1] x, bool_t[:,::1] mask):
        '''Update the running set of statistics given

        Parameters
        ----------
        x : 2d ndarray
            values from a single observation
        mask : 2d ndarray
            A uint8 array to exclude entries from the accumulated statistics.
        '''

        cdef:
            index_t i, j
            int n1
            double delta, delta_n, term1

        assert x.shape[0] == mask.shape[0] == self._sz0
        assert x.shape[1] == mask.shape[1] == self._sz1

        with nogil:
            for i in range(self._sz0):
                for j in range(self._sz1):
                    if not mask[i,j]:
                        n1 = self._n[i,j]
                        self._n[i,j] += 1
                        delta = x[i,j] - self._M1[i,j]
                        delta_n = delta / self._n[i,j]
                        term1 = delta * delta_n * n1
                        self._M1[i,j] += delta_n
                        self._M2[i,j] += term1
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __add__(StreamingStats2D self, StreamingStats2D other):
        cdef:
            index_t i, j
            int n1
            double delta, delta2
            StreamingStats2D combined 
            
        combined = StreamingStats2D((self._sz0, self._sz1))
            
        for i in range(self._sz0):
            for j in range(self._sz1):
                combined._n[i,j] = self._n[i,j] + other._n[i,j]
                delta = other._M1[i,j] - self._M1[i,j]
                delta2 = delta * delta
                combined._M1[i,j] = (other._n[i,j]*other._M1[i,j] + self._n[i,j]*self._M1[i,j]) / combined._n[i,j]
                combined._M2[i,j] = other._M2[i,j] + self._M2[i,j] + (delta2 * self._n[i,j] * other._n[i,j]) / combined._n[i,j]
        
        return combined
    
    
    def __iadd__(StreamingStats2D self, StreamingStats2D other):
        combined = self + other
        self = combined
        return self
                        

    property mean:
        def __get__(self):
            tmp = numpy.asarray(self._M1)
            return numpy.nan_to_num(tmp)

    property var:
        def __get__(self):
            tmp_m = numpy.asarray(self._M2)
            tmp_n = numpy.asarray(self._n)
            return numpy.nan_to_num(tmp_m / tmp_n)

    property n:
        def __get__(self):
            return numpy.asarray(self._n)

        def __set__(self, val):
            self._n = val[:]

    property M1:
        def __get__(self):
            return numpy.asarray(self._M1)

        def __set__(self, val):
            self._M1 = val[:]

    property M2:
        def __get__(self):
            return numpy.asarray(self._M2)

        def __set__(self, val):
            self._M2 = val[:]


cdef class StreamingStats1D:
    '''Calculate mean and variance of a series of one-dimensional arrays of shape (nbins,)
    using an online algorithm. The statistics are accumulated along what would be axis=0 if the 
    input arrays were stacked vertically. 

    This code has been adapted from:
    http://www.johndcook.com/skewness_kurtosis.html'''

    cdef weight_t[::1] _M1
    cdef weight_t[::1] _M2
    cdef uint_t[::1] _n
    cdef Py_ssize_t _sz0

    def __init__(self, int nbins):

        self._n = numpy.zeros((nbins,), dtype=numpy.uint)
        self._M1 = numpy.zeros((nbins,), dtype=weight_dtype)
        self._M2 = numpy.zeros((nbins,), dtype=weight_dtype)
        self._sz0 = nbins

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def update(self, weight_t[::1] x, bool_t[::1] mask):
        '''Update the running set of statistics given

        Parameters
        ----------
        x : 1d ndarray
            values from a single observation
        mask : 1d ndarray
            A uint8 array to exclude entries from the accumulated statistics.
        '''

        cdef:
            index_t i
            int n1
            double delta, delta_n, term1

        assert x.shape[0] == mask.shape[0] == self._sz0

        with nogil:
            for i in range(self._sz0):
                if not mask[i]:
                    n1 = self._n[i]
                    self._n[i] += 1
                    delta = x[i] - self._M1[i]
                    delta_n = delta / self._n[i]
                    term1 = delta * delta_n * n1
                    self._M1[i] += delta_n
                    self._M2[i] += term1
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __add__(StreamingStats1D self, StreamingStats1D other):
        cdef:
            index_t i
            int n1
            double delta, delta2
            StreamingStats1D combined 
            
        combined = StreamingStats1D(self._sz0)
            
        for i in range(self._sz0):
            combined._n[i] = self._n[i] + other._n[i]
            delta = other._M1[i] - self._M1[i]
            delta2 = delta * delta
            combined._M1[i] = (other._n[i]*other._M1[i] + self._n[i]*self._M1[i]) / combined._n[i]
            combined._M2[i] = other._M2[i] + self._M2[i] + (delta2 * self._n[i] * other._n[i]) / combined._n[i]
        
        return combined
    
    
    def __iadd__(StreamingStats1D self, StreamingStats1D other):
        combined = self + other
        self = combined
        return self
                        

    property mean:
        def __get__(self):
            tmp = numpy.asarray(self._M1)
            return numpy.nan_to_num(tmp)

    property var:
        def __get__(self):
            tmp_m = numpy.asarray(self._M2)
            tmp_n = numpy.asarray(self._n)
            return numpy.nan_to_num(tmp_m / tmp_n)

    property n:
        def __get__(self):
            return numpy.asarray(self._n)

        def __set__(self, val):
            self._n = val[:]

    property M1:
        def __get__(self):
            return numpy.asarray(self._M1)

        def __set__(self, val):
            self._M1 = val[:]

    property M2:
        def __get__(self):
            return numpy.asarray(self._M2)

        def __set__(self, val):
            self._M2 = val[:]
