import westpa

import logging

log = logging.getLogger(__name__)


def get_pcoord(state):
    log.debug('getting progress coordinate for {!r}'.format(state))
    propagator = westpa.rc.get_propagator()
    propagator.get_pcoord(state)
    return state


def gen_istate(basis_state, initial_state):
    log.debug('generating initial state from {!r} (into {!r})'.format(basis_state, initial_state))
    propagator = westpa.rc.get_propagator()
    propagator.update_basis_initial_states([basis_state], [initial_state])
    propagator.gen_istate(basis_state, initial_state)
    return basis_state, initial_state


def prep_iter(n_iter, segments):
    log.debug('propagator.prepare_iteration(...)')
    propagator = westpa.rc.get_propagator()
    propagator.clear_basis_initial_states()
    propagator.prepare_iteration(n_iter, segments)


def post_iter(n_iter, segments):
    log.debug('propagator.finalize_iteration(...)')
    propagator = westpa.rc.get_propagator()
    propagator.finalize_iteration(n_iter, segments)


def propagate(basis_states, initial_states, segments):
    propagator = westpa.rc.get_propagator()
    propagator.update_basis_initial_states(basis_states, initial_states)
    outgoing_ids = [segment.seg_id for segment in segments]
    incoming_segments = {segment.seg_id: segment for segment in propagator.propagate(segments)}
    if log.isEnabledFor(logging.DEBUG):
        log.debug('propagated {:d} segments'.format(len(incoming_segments)))
    return [incoming_segments[seg_id] for seg_id in outgoing_ids]
