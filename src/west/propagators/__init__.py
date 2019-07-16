# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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


import westpa
import itertools
def blocked_iter(blocksize, iterable, fillvalue = None):
    # From the Python "itertools recipes" (grouper)
    args = [iter(iterable)] * blocksize
    return itertools.zip_longest(fillvalue=fillvalue, *args)

class WESTPropagator:
    def __init__(self, rc=None):
        
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
        
        self.rc = rc or westpa.rc
        
    def prepare_iteration(self, n_iter, segments):
        """Perform any necessary per-iteration preparation.  This is run by the work manager."""
        pass
    
    def finalize_iteration(self, n_iter, segments):
        """Perform any necessary post-iteration cleanup.  This is run by the work manager."""
        pass

    # Specific functions required by the WEST framework
    def get_pcoord(self, state):
        '''Get the progress coordinate of the given basis or initial state.'''
        raise NotImplementedError
                
    def gen_istate(self, basis_state, initial_state):
        '''Generate a new initial state from the given basis state.'''
        raise NotImplementedError
                        
    def propagate(self, segments):
        """Propagate one or more segments, including any necessary per-iteration setup and teardown for this propagator."""
        raise NotImplementedError

    def clear_basis_initial_states(self):
        self.basis_states = {}
        self.initial_states = {}
    
    def update_basis_initial_states(self, basis_states, initial_states):
        self.basis_states.update({state.state_id: state for state in basis_states})
        self.initial_states.update({state.state_id: state for state in initial_states})
    