__metaclass__ = type

import itertools
def blocked_iter(blocksize, iterable, fillvalue = None):
    # From the Python "itertools recipes" (grouper)
    args = [iter(iterable)] * blocksize
    return itertools.izip_longest(fillvalue=fillvalue, *args)

class WEMDPropagator:
    def __init__(self):
        pass
        
    def prepare_iteration(self, n_iter, segments):
        """Perform any necessary per-iteration preparation.  This is run by the work manager."""
        pass
    
    def finalize_iteration(self, n_iter, segments):
        """Perform any necessary post-iteration cleanup.  This is run by the work manager."""
        pass
                        
    def propagate(self, segments):
        """Propagate one or more segments, including any necessary per-iteration setup and teardown for this propagator."""
        raise NotImplementedError
    