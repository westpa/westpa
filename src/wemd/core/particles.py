__metaclass__ = type

class Particle:
    def __repr__(self):
        return '<%s(%s): weight=%r, pcoord=%r>' \
                % (self.__class__.__name__,
                   hex(id(self)),
                   self.weight,
                   self.pcoord)
    
    def __init__(self, particle_id = None, weight = None, pcoord = None, 
                 p_parent = None, parents = None, initial_region = None):
        self.particle_id = particle_id 
        self.weight = weight
        self.pcoord = pcoord
        self.p_parent = p_parent
        self.parents = parents or set()
        self.initial_region = initial_region
        
class ParticleCollection(set):
    def __init__(self, iterable=None):
        super(ParticleCollection,self).__init__(iterable or ())
        
    def __repr__(self):
        return '<%s(%s) %d particles, norm=%g>' \
               % (self.__class__.__name__, hex(id(self)), len(self), self.norm)
        
    def get_norm(self):
        "Return the total weight of all particles in this collection"
        norm = 0.0
        for particle in self:
            norm += particle.weight
        return norm
            
    norm = property(get_norm, None, None, 
                    "The total weight of all particles in this collection")

    def renorm(self, weight):
        current_weight = self.get_norm()
        if current_weight == 0.0 and weight != 0.0:
            raise ValueError('cannot renormalize empty bin')
        else:
            wfactor = weight / current_weight
            for particle in self:
                particle.weight *= wfactor
