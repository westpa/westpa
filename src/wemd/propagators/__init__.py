__metaclass__ = type

import itertools
def blocked_iter(blocksize, iterable, fillvalue = None):
    # From the Python "itertools recipes" (grouper)
    args = [iter(iterable)] * blocksize
    return itertools.izip_longest(fillvalue=fillvalue, *args)

class WEMDPropagator:
    def __init__(self):
        
        # For maximum flexibility, the basis states and initial states valid
        # at the point in the simulation when the propgator is used must be
        # available in several routines, and it is inconvenient to pass them
        # to every routine that needs them. A currently-reasonable-seeming solution
        # is to store at least the basis states and initial states necessary for
        # the current operation (propagation, etc). The set_basis_initial_states() function 
        # accomplishes this. They are stored as dictionaries of state_id -> state,
        # so they can be looked up by ID without needing to store them all (and 
        # thus potentially send them all over the wire when only one of them is needed, e.g.)        
        self.basis_states = {}
        self.initial_states = {}
        
    def prepare_iteration(self, n_iter, segments):
        """Perform any necessary per-iteration preparation.  This is run by the work manager."""
        pass
    
    def finalize_iteration(self, n_iter, segments):
        """Perform any necessary post-iteration cleanup.  This is run by the work manager."""
        pass
                        
    def propagate(self, segments):
        """Propagate one or more segments, including any necessary per-iteration setup and teardown for this propagator."""
        raise NotImplementedError

    def clear_basis_initial_states(self):
        self.basis_states = {}
        self.initial_states = {}
    
    def update_basis_initial_states(self, basis_states, initial_states):
        self.basis_states.update({state.state_id: state for state in basis_states})
        self.initial_states.update({state.state_id: state for state in initial_states})
    