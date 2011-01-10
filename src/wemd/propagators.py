__metaclass__ = type

import itertools
def blocked_iter(blocksize, iterable, fillvalue = None):
    # From the Python "itertools recipes" (grouper)
    args = [iter(iterable)] * blocksize
    return itertools.izip_longest(fillvalue=fillvalue, *args)


class WEMDPropagator:
    def __init__(self, work_manager):
        # Reference to the parent work manager
        self.work_manager = work_manager
                
    def pre_iter(self, n_iter):
        pass
    
    def post_iter(self, n_iter):
        pass
    
    def pre_segment(self, segment):
        pass
    
    def post_segment(self, segment):
        pass
    
    def _propagate_segment(self, segment):
        raise NotImplementedError
    
    def propagate_segments(self, segments):
        for segment in segments:
            self.pre_segment(segment)
            self._propagate_segment(segment)
            self.post_segment(segment)
        
     
    
    