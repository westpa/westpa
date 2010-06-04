import os, sys, subprocess, time, datetime, tempfile, string, re
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
                                    'backend.test.nsteps',
                                    'bins.ndim',
                                    'wemd.initial_pcoord'
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
        self._dimension = sim_config['bins.ndim'] = sim_config_src.get_int('bins.ndim')
        
        # number of timesteps
        self._nsteps = sim_config['backend.test.nsteps'] = sim_config_src.get_int('backend.test.nsteps')
        
        #read in intial coord, in nm
        if 'wemd.initial_pcoord' in sim_config.keys():
            self._initial_pcoord = sim_config['wemd.initial_pcoord']
        else:
            pcoord_vals = [float(x) for x in sim_config_src.get_list('wemd.initial_pcoord')]
            self._initial_pcoord = sim_config['wemd.initial_pcoord'] = numpy.array( pcoord_vals, dtype=numpy.float64 )
        
        #read in the Gaussian Potentials (if present)
        potential = {}
        potential_entries = [key for key in sim_config_src if key.startswith('backend.test.potential')]
        reIsPotential = re.compile('backend\.test\.potential_(.+)_(.+)')
        for potential_entry in potential_entries:
            m = reIsPotential.match( potential_entry )
            if not m:
                raise ConfigError("Invalid potential specified")
            else:
                if m.group(2) not in potential.keys():
                    potential[m.group(2)] = {}
                    
                if m.group(1) == 'height':
                    potential[m.group(2)][m.group(1)] = sim_config_src.get_float(potential_entry)
                elif m.group(1) == 'coord' or m.group(1) == 'width':
                    coord_vals = [float(x) for x in sim_config_src.get_list(potential_entry)]
                    potential[m.group(2)][m.group(1)] = numpy.array( coord_vals, dtype=numpy.float64 )
                else:
                    raise ConfigError('Invalid potential specified')

        #check to make sure they have the correct dimensions
        for pkey in potential.keys():
            for key in ('height', 'coord', 'width'):
                if key not in potential[pkey]:
                    raise ConfigError("Invalid potential -- missing key %s" % (key))
            
            if len(potential[pkey]['coord']) != self._dimension or len(potential[pkey]['width']) != self._dimension:
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
            for key in ('temperature', 'timestep', 'mass', 'friction_const',
                        'potentials', 'nsteps'):
                setattr(self, '_%s' % key, self.sim_config['backend.test.%s' % key])
                
            self._initial_pcoord = self.sim_config['wemd.initial_pcoord']
            self._dimension = self.sim_config['bins.ndim']
            self._source_pcoords = self.sim_config['bin.source_pcoords']

            
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            self.pre_segment(segment)
                       
            if segment.n_iter == 1 or segment.p_parent == None:
                log.debug('Initializing segment coordinates from initial pcoord')

                if segment.data.get('initial_region') is not None:
                    initial_region = segment.data['initial_region']

                    source_pcoord = self._source_pcoords[initial_region]
                    source_pcoord_val = numpy.array( [source_pcoord['pcoord']], dtype=numpy.float64 )
                    
                    #use the pcoord associated with the region
                    segment.pcoord = copy(source_pcoord_val)                    

                else:
                    pcoord_val = numpy.array( [self._initial_pcoord], dtype=numpy.float64 )
                    segment.pcoord = copy(pcoord_val)
           
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

        for pkey in self._potentials.keys():       
    
            # H (k_b T) * nm
            g_const = ( 1 / nano ) * self._potentials[pkey]['height'] * k_b * self._temperature

            # exp(- (x-mu_x)^2/(2 sigma_x^2)) * ...
            g_exp = 1
            for i in xrange(0, dim):
                g_exp *= math.exp(-0.5 * (pos[i] - self._potentials[pkey]['coord'][i]) ** 2.0 / (self._potentials[pkey]['width'][i] ** 2.0) )

            g_grad = []
            for i in xrange(0, dim):
                g_grad.append( - g_const * g_exp * (pos[i] - self._potentials[pkey]['coord'][i]) / (self._potentials[pkey]['width'][i] ** 2.0) )

            diff_potential_val += numpy.array(g_grad,dtype=numpy.float64)

        return diff_potential_val
