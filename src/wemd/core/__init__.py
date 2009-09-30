class ConfigError(ValueError):
    def __init__(self, message = '', exc = None):
        self.message = message
        self.exc = exc

    def __str__(self):
        if self.message:
            if self.exc:
                return '%s: %s' % (self.message, self.exc)
            else:
                return self.message
        elif self.exc:
            return str(self.exc)
        else:
            return ValueError.__str__(self)

import particles, segments, we_sim
from particles import Particle, ParticleCollection
from segments import Segment
from we_sim import WEError, ParticleExitError