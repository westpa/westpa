'''
Created on Feb 11, 2012

@author: mzwier
'''

import logging
log = logging.getLogger(__name__)

def get_pcoord(propagator, state):
    propagator.get_pcoord(state)
    return state
    
def gen_istate(propagator, basis_state, initial_state):
    propagator.gen_istate(basis_state, initial_state)
    return basis_state, initial_state
    
def prep_iter(propagator, n_iter, segments):
    propagator.prepare_iteration(n_iter, segments)
    
def post_iter(propagator, n_iter, segments):
    propagator.prepare_iteration(n_iter, segments)
    
def propagate(propagator, segments):
    '''Propagate the given segments with the given propagator. This has to be a top-level
    function for the current incarnation of the work manager.'''
    outgoing_ids = [segment.seg_id for segment in segments]
    incoming_segments = {segment.seg_id: segment for segment in propagator.propagate(segments)}
    if log.isEnabledFor(logging.DEBUG):
        log.debug('propagated {!r}'.format(incoming_segments))
    return [incoming_segments[seg_id] for seg_id in outgoing_ids]
