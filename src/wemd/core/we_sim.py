__metaclass__ = type

class WESimIter:
    """
    Describes per-iteration information (summary or otherwise) for
    a WE simulation.
    """
    
    def __init__(self, i_iter = None, n_particles = None, norm = None,
                 cputime = None, walltime = None, binarray = None,
                 segments = None):
        self.i_iter = i_iter
        self.n_particles = n_particles
        self.norm = norm
        self.cputime = cputime
        self.walltime = walltime
        self.binarray = binarray