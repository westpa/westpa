import os, sys, subprocess, time, datetime, tempfile, string
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
from scipy.spatial.distance import pdist
import cPickle

from wemd.core import Segment
from wemd.util.mpi import getrank
from wemd.backend_drivers import BackendDriver
from DoubleGo import DoubleGo

import simtk
import simtk.chem.openmm as openmm
import simtk.unit as units

class PyOpenMMBackend(BackendDriver):
    """
    Backend driver to run Langevin Dynamics using PyOpenMM
    All parameters and setting are passed from the runtime configuration file and
    must include the specification of a class that takes a single configuration file
    and has methods: 
    getSystem() - returns a valid openmm system
    getCoordinates(int state) where state is 0 or 1 for the initial and target states
        respectively. This method must return valid numpy arrays of type units.Quantity
        
    By default this driver uses the dRMSD between the current and target state
    as the progress coordinate.
    """
    
    def __init__(self, runtime_config):
        super(PyOpenMMBackend, self).__init__(runtime_config)
        
        assert(runtime_config['backend.driver'].lower() == 'pyopenmm')
        
        # Get the system creator object and configuration file
        sysmaker = eval(runtime_config.get('backend.pyopenmm.sysmaker'))
        sysconfig = runtime_config.get('backend.pyopenmm.sysconfig')
        
        # Setup and extract the system
        self._OpMM = sysmaker(sysconfig)
        self._system = self._OpMM.getSystem()
        
        self._steps_per_seg = runtime_config.get_int('backend.pyopenmm.steps_per_segment')
        log.debug("Will run %d steps per segment" % (self._steps_per_seg))
        
        self._segdir = runtime_config.get('backend.pyopenmm.segdir')
        
        pp = runtime_config.get('backend.pyopenmm.platform')
        self._platform = openmm.Platform.getPlatformByName(pp)
        
        # Define Langevin Integrator parameters
        self._temperature = eval(runtime_config.get('backend.pyopenmm.temperature'))
        self._collision_rate = eval(runtime_config.get('backend.pyopenmm.collision_rate'))
        self._timestep = eval(runtime_config.get('backend.pyopenmm.timestep'))
        
        # Get initial and target coordinates
        self._initial_coord = self._OpMM.getCoordinates(0)
        self._target_coord = self._OpMM.getCoordinates(1)
        self._initial_pdist = pdist(self._initial_coord,'euclidean')
        self._target_pdist = pdist(self._target_coord,'euclidean')
        N = self._target_pdist.shape[0]
        self._pcfactor = 1.0/N
        
        ipcoord = -numpy.sqrt(self._pcfactor*numpy.sum((self._target_pdist - self._initial_pdist)**2))
        log.debug("Initial pcoord = %f" % (ipcoord))
        
        
    def pre_iter(self, we_iter):
        pass
    
    def post_iter(self, we_iter):
        pass
    
    def pre_segment(self, segment):
        pass
    
    def post_segment(self, segment):
        pass
    
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            self.pre_segment(segment)

            integrator = openmm.LangevinIntegrator(self._temperature, self._collision_rate, self._timestep)
            rank = getrank()
            rn = int(time.time()) / (rank + 1) + rank
            integrator.setRandomNumberSeed(rn)
            context = openmm.Context(self._system, integrator, self._platform)
            
            if segment.n_iter == 1 or segment.p_parent == None:
                log.debug('Initializing segment coordinates from reference structure')
                coords = self._OpMM.getCoordinates(0)
                context.setPositions(coords)
            else:
                log.debug('Initializing segment coordinates from parent structure')
                coordpkl = self._segdir + "/%d/%d/" % (segment.p_parent.n_iter, segment.p_parent.seg_id) + "coord.pkl"
                cpklfile = open(coordpkl,'rb')
                coords = cPickle.load(cpklfile)
                cpklfile.close()
                
                velpkl = self._segdir + "/%d/%d/" % (segment.p_parent.n_iter, segment.p_parent.seg_id) + "vel.pkl"
                vpklfile = open(velpkl,'rb')
                velocs = cPickle.load(vpklfile)
                vpklfile.close()
                
                context.setPositions(coords)
                context.setVelocities(velocs)
                      
            # Record start timing info
            segment.starttime = datetime.datetime.now()
            init_walltime = time.time()
            init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            log.debug('launched at %s' % segment.starttime)

            # Propagate System
            integrator.step(self._steps_per_seg)

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

            state = context.getState(getPositions=True,getVelocities=True)
            newcoords = state.getPositions() 
            newvelocs = state.getVelocities() 

            self._update_progress_coordinate(segment,newcoords)
            
            # Dump coordinates and velocities into cPickle
            bdir = self._segdir + "/%d/%d/" % (segment.n_iter, segment.seg_id)
            os.makedirs(bdir)
            coordpkl = bdir + "/coord.pkl"
            velpkl = bdir + "/vel.pkl"
            
            cpklfile = open(coordpkl,'wb')
            vpklfile = open(velpkl,'wb')
            
            cPickle.dump(newcoords / units.nanometer,cpklfile,cPickle.HIGHEST_PROTOCOL)
            cPickle.dump(newvelocs / (units.nanometer/units.picosecond),vpklfile,cPickle.HIGHEST_PROTOCOL)
            
            cpklfile.close()
            vpklfile.close()
            
            # Clean up
            del integrator
            del context
            del state   

            log.debug('Segment run successful')
            segment.status = Segment.SEG_STATUS_COMPLETE
                            
            self.post_segment(segment)

    def _update_progress_coordinate(self,segment,coord):
        """Calculate dRMSD from coordinates"""
        npd = pdist(coord.value_in_unit(units.angstrom),'euclidean')
        pcoord = -numpy.sqrt(self._pcfactor*numpy.sum((self._target_pdist - npd)**2))
        log.debug('pcoord = %f' % (pcoord))
        segment.pcoord = numpy.array([[pcoord]], dtype=numpy.float64)



         
     


