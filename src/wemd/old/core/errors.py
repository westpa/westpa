from wemd.util.config_dict import ConfigError

class WEError(Exception):
    pass

class WEEnvironmentError(WEError, EnvironmentError):
    pass

class WEEnvironmentVarError(WEEnvironmentError):
    def __init__(self, e):
        self.message = 'WE environment variable not set: %s' % e
        
class PropagationIncompleteError(WEError):
    pass

class ParticleExitError(WEError):
    """Particle has exited the bin space
    args[0] = particle
    args[1] = tuple index into the bin limit array; None appears where the 
              particle did not fall into specified limits
              
    This may be raised as an error condition, or as a (rather inefficient)
    signal to process the particle in some way (after which the swarm will
    have to be re-binned from scratch)."""
    
    def __init__(self, particle, index):
        self.particle = particle
        self.index = index

    def __str__(self):
        return 'particle exited bin space; pcoord=%s, indices=%s' \
               % (self.particle.pcoord, self.index)  

class WEConfigError(ConfigError): pass
