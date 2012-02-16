'''
Created on Feb 11, 2012

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

def get_pcoord(propagator, state):
    log.debug('getting progress coordinate for {!r}'.format(state))
    propagator.get_pcoord(state)
    return state
    
def gen_istate(propagator, basis_state, initial_state):
    log.debug('generating initial state from {!r} (into {!r})'.format(basis_state, initial_state))
    propagator.gen_istate(basis_state, initial_state)
    return basis_state, initial_state
    
def prep_iter(propagator, n_iter, segments):
    log.debug('propagator.prepare_iteration(...)')
    propagator.prepare_iteration(n_iter, segments)
    
def post_iter(propagator, n_iter, segments):
    log.debug('propagator.finalize_iteration(...)')
    propagator.finalize_iteration(n_iter, segments)
    
def propagate(propagator, basis_states, initial_states, segments):
    '''Propagate the given segments with the given propagator. This has to be a top-level
    function for the current incarnation of the work manager.'''
    propagator.update_basis_initial_states(basis_states, initial_states)
    outgoing_ids = [segment.seg_id for segment in segments]
    incoming_segments = {segment.seg_id: segment for segment in propagator.propagate(segments)}
    if log.isEnabledFor(logging.DEBUG):
        log.debug('propagated {!r}'.format(incoming_segments))
    return [incoming_segments[seg_id] for seg_id in outgoing_ids]
