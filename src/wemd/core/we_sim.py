__metaclass__ = type

class WESimIter:
    """
    Describes per-iteration information (summary or otherwise) for
    a WE simulation.
    """
    
    def __init__(self, we_iter = None, n_particles = None, norm = None,
                 cputime = None, walltime = None, binarray = None):
        self.we_iter = we_iter
        self.n_particles = n_particles
        self.norm = norm
        self.cputime = cputime
        self.walltime = walltime
        self.binarray = binarray