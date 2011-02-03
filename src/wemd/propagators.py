__metaclass__ = type

import itertools
def blocked_iter(blocksize, iterable, fillvalue = None):
    # From the Python "itertools recipes" (grouper)
    args = [iter(iterable)] * blocksize
    return itertools.izip_longest(fillvalue=fillvalue, *args)

class WEMDPropagator:
    def __init__(self, sim_manager):
        # Reference to the parent work manager
        self.sim_manager = sim_manager
        
    def prepare_iteration(self, n_iter):
        pass
    
    def finalize_iteration(self, n_iter):
        pass
                    
    def prepare_segment(self, segment):
        pass
            
    def finalize_segment(self, segment):
        pass
    
    def propagate_segments(self, segments):
        raise NotImplementedError
    