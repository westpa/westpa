import os, sys, subprocess, time, datetime, tempfile, string
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
from wemd.core import Segment
from wemd.util.config_dict import ConfigError
from wemd.backend_drivers import BackendDriver

import math
from copy import copy

k_b = 1.3806503 * 10 ** -23.0 # boltzmann const m^2 kg s^-2 K^-1
nano = 10.0 ** -9.0
class TestBackend(BackendDriver):
    """
    Backend driver to run a test system using discrete overdamped Langevin Dynamics
    """

    def __init__(self):
        super(TestBackend, self).__init__()

        
    def sim_init(self, sim_config, sim_config_src):
        self.sim_config = sim_config
        sim_config_src.require_all(('backend.test.temperature',
                                    'backend.test.timestep',
                                    'backend.test.mass',
                                    'backend.test.friction_const',
                                    'backend.test.dimension',
                                    'backend.test.initial_pcoord',
                                    'backend.test.nsteps'
                                    ))
        # temperature in K
        self._temperature = sim_config['backend.test.temperature'] = sim_config_src.get_float('backend.test.temperature')
        #timestep in ps
        self._timestep = sim_config['backend.test.timestep'] = sim_config_src.get_float('backend.test.timestep') * 10.0 ** -12.0
        # mass in kg
        self._mass = sim_config['backend.test.mass'] = sim_config_src.get_float('backend.test.mass')
        # friction constant (1/s)
        self._friction_const = sim_config['backend.test.friction_const'] = sim_config_src.get_float('backend.test.friction_const')
        # dimensionality
        self._dimension = sim_config['backend.test.dimension'] = sim_config_src.get_int('backend.test.dimension')
        # number of timesteps
        self._nsteps = sim_config['backend.test.nsteps'] = sim_config_src.get_int('backend.test.nsteps')
        
        #read in intial coord, in nm
        pcoord_vals = sim_config_src.get_list('backend.test.initial_pcoord', type=float)
        self._initial_pcoord = sim_config['backend.test.initial_pcoord'] = numpy.array( [pcoord_vals], dtype=numpy.float64 )

        #read in the Gaussian Potentials (if present)
        potential = []
        potential_entries = [key for key in sim_config_src if key.startswith('backend.test.potential')]
        for potential_entry in potential_entries:
            # FIXME: prefer typed retrieval functions from ConfigDict over eval
            potential.append(eval(sim_config_src.get(potential_entry)))

        #check to make sure they have the correct dimensions
        for i in xrange(0, len(potential)):
            if len(potential[i][1]) != self._dimension or len(potential[i][2]) != self._dimension:
                raise ConfigError("Invalid potential specified, check dimensions")

        self._potentials = sim_config['backend.test.potentials'] = potential

        #allow for > 3 dimensions
        if( self._dimension < 1 ):
            raise ConfigError("Invalid Dimension")
        
    def runtime_init(self, runtime_config):
        super(TestBackend, self).runtime_init(runtime_config)
        
        try:
            self.sim_config['backend.test.temperature']
        except (KeyError,TypeError):
            log.info('skipping backend configuration for now')
        else:            
            for key in ('temperature', 'timestep', 'mass', 'friction_const', 'dimension',
                        'initial_pcoord', 'potentials', 'nsteps'):
                setattr(self, '_%s' % key, self.sim_config['backend.test.%s' % key])
    
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            self.pre_segment(segment)

        
            if segment.n_iter == 1 or segment.p_parent == None:
                log.debug('Initializing segment coordinates from initial pcoord')
                
                if segment.data.get('initial_region') is not None:
                    segment.pcoord = copy(self._initial_pcoord)                    
                    print "test backend using region %r" % segment.data['initial_region']
                else:
                    segment.pcoord = copy(self._initial_pcoord)
            else:
                log.debug('Initializing segment coordinates from parent')
                segment.pcoord = copy(segment.p_parent.pcoord)

            # Record start timing info
            segment.starttime = datetime.datetime.now()
            init_walltime = time.time()
            init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            log.debug('launched at %s' % segment.starttime)

            segment.data['dt'] = self._timestep
            # Propagate System
            self.update_progress_coordinate(segment)

            #print "segment.pcoord: %r" % segment.pcoord

            # Record end timing info
            final_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            final_walltime = time.time()
            segment.endtime = datetime.datetime.now()
            segment.walltime = final_walltime - init_walltime
            segment.cputime = final_cputime - init_cputime
            log.debug('completed at %s (wallclock %s, cpu %s)'
                      % (segment.endtime,
                         segment.walltime,
                         segment.cputime))

            log.debug('Segment run successful')
            segment.status = Segment.SEG_STATUS_COMPLETE
                            
            self.post_segment(segment)

    def update_progress_coordinate(self,segment):
        
        # tau / ( m * gamma ) (gamma is friction constant)
        tmg = self._timestep / ( self._mass * self._friction_const )

        #variance for Gaussian
        var = 2.0 * k_b * self._temperature * tmg
        std = math.sqrt( var )
        dim = self._dimension
        nsteps = self._nsteps

        #extend pcoord by the number of integration steps
        pcoord = numpy.empty((segment.pcoord.shape[0]+nsteps,segment.pcoord.shape[1]),dtype=numpy.float64)
      
        #copy the existing pcoords
        pcoord_size = len(segment.pcoord)
        for i in xrange(0, pcoord_size):
            pcoord[i] = segment.pcoord[i]

        for itime in xrange(pcoord_size, pcoord_size+nsteps):
            #discrete overdamped Langevin equation
            #x_{j+1} = x_j - {delta t} over {m %gamma} (dU/dx)_xj + delta{x^rand}
            pcoord[itime] = pcoord[itime-1] - tmg * self.diff_potential( pcoord[itime-1], dim ) / nano + numpy.random.normal(0, std, (1,dim)) / nano

        segment.pcoord = copy(pcoord)

    def diff_potential(self,pos,dim):
        """returns the gradient of the potential as a numpy.array
           must give array appropriate for dimension
           in units of N
        """

        #no potential present
        if not self._potentials:
            return numpy.array([0 for i in range(0,dim)],dtype=numpy.float64)

        diff_potential_val = numpy.zeros((1,dim),dtype=numpy.float64)

        #potential has the format [ H (kT), [mu_x,...] (nm), [sigma_x,...] (nm) ]
        for potential in self._potentials:

            # H (k_b T) * nm
            g_const = ( 1 / nano ) * potential[0] * k_b * self._temperature

            # exp(- (x-mu_x)^2/(2 sigma_x^2)) * ...
            g_exp = 1
            for i in xrange(0, dim):
                g_exp *= math.exp(-0.5 * (pos[i] - potential[1][i]) ** 2.0 / (potential[2][i] ** 2.0) )

            g_grad = []
            for i in xrange(0, dim):
                g_grad.append( - g_const * g_exp * (pos[i] - potential[1][i]) / (potential[2][i] ** 2.0) )

            diff_potential_val += numpy.array(g_grad,dtype=numpy.float64)

        return diff_potential_val
