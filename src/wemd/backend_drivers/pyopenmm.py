import os, sys, subprocess, time, datetime, tempfile, string
import errno,exceptions
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
from scipy.spatial.distance import pdist
import cPickle

from wemd.core import Segment
from wemd.util.mpi import getrank
from wemd.backend_drivers import BackendDriver

import simtk
import simtk.chem.openmm as openmm
import simtk.unit as units

#=============================================================================================
# MODULE CONSTANTS
#=============================================================================================
kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant
#=============================================================================================

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
        
        #assert(runtime_config['backend.driver'].lower() == 'pyopenmm')
        
        # Get the system creator object and configuration file
        system_module_name = runtime_config.get('backend.pyopenmm.system_module')
        
        from wemd.core import ConfigError
        try:
            system_module = sys.modules[system_module_name] = __import__(system_module_name,globals(),locals(),[],-1)
        except ImportError, e:
            raise ConfigError('pyopenmm backend driver system module unavailable (%s)' % e)
        
        try:
            sysmaker = getattr(system_module,runtime_config.get('backend.pyopenmm.sysmaker'))
        except AttributeError, e:
            raise ConfigError('pyopenmm backend driver system class unavailable (%s)' % e)
        
        
        sysconfig = runtime_config.get('backend.pyopenmm.sysconfig')
        
        # Setup and extract the system
        self._OpMM = sysmaker(sysconfig)
        self._system = self._OpMM.getSystem()
        
        self._steps_per_seg = runtime_config.get_int('backend.pyopenmm.steps_per_segment')
        log.debug("Will run %d steps per segment" % (self._steps_per_seg))
        
        self._segdir = runtime_config.get('backend.pyopenmm.segdir')
        
        pp = runtime_config.get('backend.pyopenmm.platform')
        self._platform = openmm.Platform.getPlatformByName(pp)
        
        # TODO: Need to devise a better way of handling device assignment when using multiple hosts
        if pp == 'OpenCL':
            self._platform.setPropertyDefaultValue("OpenCLDeviceIndex",str(getrank()))
        elif pp == 'Cuda':
            self._platform.setPropertyDefaultValue("CudaDevice",str(getrank()))
            
        # Define Langevin Integrator parameters
        self._temperature = eval(runtime_config.get('backend.pyopenmm.temperature'))
        self._collision_rate = eval(runtime_config.get('backend.pyopenmm.collision_rate'))
        self._timestep = eval(runtime_config.get('backend.pyopenmm.timestep'))
        
        self._numParticles = self._OpMM.getNumParticles()
        
        # Setup context and integrator
        self._integrator = openmm.LangevinIntegrator(self._temperature, self._collision_rate, self._timestep)
        rank = getrank()
        rn = int(time.time()) / (rank + 1) + rank
        self._integrator.setRandomNumberSeed(rn)
        self._context = openmm.Context(self._system, self._integrator, self._platform)
        
        # Compute thermal energy and inverse temperature from specified temperature.
        self._kT = kB * self._temperature # thermal energy
        self._beta = 1.0 / self._kT # inverse temperature
        self._mass = units.Quantity(numpy.zeros([self._numParticles, 1], numpy.float32), units.amu)
        
        self._sigma = units.Quantity(numpy.zeros([self._numParticles, 1], numpy.float32), units.sqrt(self._kT.unit/units.amu))
        for atom_index in xrange(self._numParticles):
            self._mass[atom_index] = self._system.getParticleMass(atom_index)
            self._sigma[atom_index] = units.sqrt(self._kT/self._mass[atom_index])
        
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
            
            if segment.n_iter == 1 or segment.p_parent == None:
                log.debug('Initializing segment coordinates from reference structure')
                coords = self._OpMM.getCoordinates(0)
                velocs = units.Quantity(numpy.random.normal(size=(self._numParticles,3)) * self._sigma,self._sigma.unit)
                
                self._context.setPositions(coords)
                self._context.setVelocities(velocs)
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
                
                self._context.setPositions(coords)
                self._context.setVelocities(velocs)
                      
            # Record start timing info
            segment.starttime = datetime.datetime.now()
            init_walltime = time.time()
            init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            log.debug('launched at %s' % segment.starttime)

            # Propagate System
            self._integrator.step(self._steps_per_seg)

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

            state = self._context.getState(getPositions=True,getVelocities=True)
            newcoords = state.getPositions() 
            newvelocs = state.getVelocities() 

            self._update_progress_coordinate(segment,newcoords)
            
            # Dump coordinates and velocities into cPickle
            bdir = self._segdir + "/%d/%d/" % (segment.n_iter, segment.seg_id)
            self.make_dirs(bdir)
            coordpkl = bdir + "/coord.pkl"
            velpkl = bdir + "/vel.pkl"
            
            cpklfile = open(coordpkl,'wb')
            vpklfile = open(velpkl,'wb')
            
            cPickle.dump(newcoords / units.nanometer,cpklfile,cPickle.HIGHEST_PROTOCOL)
            cPickle.dump(newvelocs / (units.nanometer/units.picosecond),vpklfile,cPickle.HIGHEST_PROTOCOL)
            
            cpklfile.close()
            vpklfile.close() 

            log.debug('Segment run successful')
            segment.status = Segment.SEG_STATUS_COMPLETE
                            
            self.post_segment(segment)

    def _update_progress_coordinate(self,segment,coord):
        """Calculate dRMSD from coordinates"""
        npd = pdist(coord.value_in_unit(units.angstrom),'euclidean')
        pcoord = -numpy.sqrt(self._pcfactor*numpy.sum((self._target_pdist - npd)**2))
        log.debug('pcoord = %f' % (pcoord))
        segment.pcoord = numpy.array([[pcoord]], dtype=numpy.float64)
        

    def make_dirs(self,dirname):
        tx = None
        try:
            os.makedirs(dirname)
        except OSError,x:
            tx = x

        if not os.path.isdir(dirname):
            if tx:
                raise tx
            raise exceptions.IOError, "unknown error prevented the creation of directory: %s" % dirname
        
class PyOpenMMBackendMultiSeg(PyOpenMMBackend):
    """docstring for PyOpenMMBackendMultiSeg"""
    
    def __init__(self, runtime_config):
        super(PyOpenMMBackendMultiSeg, self).__init__(runtime_config)
        
        # Check to make sure that the block size is equal to the number of replicas
        self._blockSize = runtime_config.get_int('backend.blocksize', 1)
        if self._blockSize != self._OpMM.getNumReplicas():
            raise AssertionError("The segment block size (%d) does not match the number of replicas in the system (%d)" \
                                    % (self._blockSize,self._OpMM.getNumReplicas()))
                                    
        self._grid = self._OpMM.getGrid(self._blockSize)
        self._gridSpacing = eval(runtime_config.get('backend.pyopenmm.gridspacing'))
        (self._xL,self._yL,self._zL) = self._OpMM.getDimensions()
                          
        
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        if len(segments) != self._OpMM.getNumReplicas():
            raise AssertionError("The number of segments (%d) does not match the number of replicas in the system (%d)" \
                                    % (len(segments),self._OpMM.getNumReplicas()))
                                    
        if not all(s != None for s in segments):
            raise AssertionError("Attempted to pass incomplete workload to PyOpenMMMultiSeg propagate_segments")
        
        coords = None
        velocs = None
        replicaCount = 0

        for segment in segments:
            # For each segment determine if initializing from parent or reference
            if segment.n_iter == 1 or segment.p_parent == None:
                log.debug('Initializing segment coordinates from reference structure')
                segcoords = self._OpMM.getCoordinates(0)[0:self._numParticles,:] / units.nanometer
                segvelocs = units.Quantity(numpy.random.normal(size=(self._numParticles,3)) * self._sigma,self._sigma.unit)
                segvelocs = segvelocs / (units.nanometer/units.picosecond)
                
            else:
                log.debug('Initializing segment coordinates from parent structure')
                coordpkl = self._segdir + "/%d/%d/" % (segment.p_parent.n_iter, segment.p_parent.seg_id) + "coord.pkl"
                cpklfile = open(coordpkl,'rb')
                segcoords = cPickle.load(cpklfile)
                cpklfile.close()
                
                velpkl = self._segdir + "/%d/%d/" % (segment.p_parent.n_iter, segment.p_parent.seg_id) + "vel.pkl"
                vpklfile = open(velpkl,'rb')
                segvelocs = cPickle.load(vpklfile)
                vpklfile.close()
            
            (NX,NY,NZ) = self._OpMM.getGridIndex(replicaCount)
        
            # Find the center-of-geometry of the seg's coords and translate to origin
            CoG = segcoords.sum(axis=0) / self._numParticles
            segcoords = segcoords - numpy.tile(CoG,(self._numParticles,1))
            
            # Move coords to lattice position
            segcoords[:,0] = segcoords[:,0] + (NX-1)*((self._xL + self._gridSpacing) / units.nanometer)
            segcoords[:,1] = segcoords[:,1] + (NY-1)*((self._yL + self._gridSpacing) / units.nanometer)
            segcoords[:,2] = segcoords[:,2] + (NZ-1)*((self._zL + self._gridSpacing) / units.nanometer) 
            
            # Concatemate coordinates to lattice coords
            if coords == None:    
                coords = segcoords.copy()
                velocs = segvelocs.copy()
            else:
                coords = units.Quantity(numpy.concatenate((coords,segcoords)), units.nanometer)
                velocs = units.Quantity(numpy.concatenate((velocs,segvelocs)), units.nanometer/units.picosecond)
        
            replicaCount += 1
        
        self._context.setPositions(coords)
        self._context.setVelocities(velocs)
            
        # Record start timing info
        for segment in segments:
            segment.starttime = datetime.datetime.now()
            log.debug('launched at %s' % segment.starttime)
        init_walltime = time.time()
        init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
        

        # Propagate System
        self._integrator.step(self._steps_per_seg)

        # Record end timing info
        final_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
        final_walltime = time.time()
        for segment in segments:
            segment.endtime = datetime.datetime.now()
            segment.walltime = final_walltime - init_walltime
            segment.cputime = final_cputime - init_cputime
            log.debug('completed at %s (wallclock %s, cpu %s)'
                  % (segment.endtime,
                     segment.walltime,
                     segment.cputime))

        state = self._context.getState(getPositions=True,getVelocities=True)
        newcoords = state.getPositions(asNumpy=True) 
        newvelocs = state.getVelocities(asNumpy=True) 

        replicaCount = 0

        for segment in segments:
            segcoords = newcoords[replicaCount*self._numParticles:(replicaCount+1)*self._numParticles,:]
            segvelocs = newvelocs[replicaCount*self._numParticles:(replicaCount+1)*self._numParticles,:]

            # Determine progress coordinate for each replica 
            self._update_progress_coordinate(segment,segcoords)
            
            # Dump each segment's new position and velocity into pickle
            bdir = self._segdir + "/%d/%d/" % (segment.n_iter, segment.seg_id)
            self.make_dirs(bdir)
            coordpkl = bdir + "/coord.pkl"
            velpkl = bdir + "/vel.pkl"
            
            cpklfile = open(coordpkl,'wb')
            vpklfile = open(velpkl,'wb')
            
            cPickle.dump(segcoords / units.nanometer,cpklfile,cPickle.HIGHEST_PROTOCOL)
            cPickle.dump(segvelocs / (units.nanometer/units.picosecond),vpklfile,cPickle.HIGHEST_PROTOCOL)
            
            cpklfile.close()
            vpklfile.close()
        
            segment.status = Segment.SEG_STATUS_COMPLETE
        
            replicaCount += 1 

        log.debug('Segments run successful')
        
             
            
     

     
     


