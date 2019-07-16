import os
import errno
import random
import time
import numpy as np
import west
from west.propagators import WESTPropagator
from west import Segment
from west.states import BasisState, InitialState

import simtk.openmm.openmm as openmm
import simtk.unit as units

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

pcoord_len = 11
pcoord_dtype = np.float32


class OpenMMPropagator(WESTPropagator):
    def __init__(self, rc=None):
        super(OpenMMPropagator, self).__init__(rc)

        self.pcoord_len = pcoord_len
        self.pcoord_dtype = pcoord_dtype
        self.pcoord_ndim = 1

        self.basis_coordinates = np.array([[5.0, 0.0, 0.0], [-5.0, 0.0, 0.0]], dtype=pcoord_dtype)

        # Default platform properties
        self.platform_properties = {'OpenCLPrecision': 'mixed',
                                   'OpenCLPlatformIndex': '0',
                                   'OpenCLDeviceIndex': '0',
                                   'CudaPrecision': 'mixed',
                                   'CudaDeviceIndex': '0'}

        config = self.rc.config

        # Validate configuration
        for key in [('west', 'openmm', 'system', 'file'),
                    ('west', 'openmm', 'integrator', 'file'),
                    ('west', 'openmm', 'integrator', 'steps_per_tau'),
                    ('west', 'openmm', 'integrator', 'steps_per_write'),
                    ('west', 'openmm', 'platform', 'name'),
                    ('west', 'data', 'data_refs', 'initial_state')]:
            config.require(key)

        self.initial_state_ref_template = config['west','data','data_refs','initial_state']

        system_xml_file = config['west', 'openmm', 'system', 'file']
        self.integrator_xml_file = config['west', 'openmm', 'integrator', 'file']

        self.steps_per_tau = config['west', 'openmm', 'integrator', 'steps_per_tau']
        self.steps_per_write = config['west', 'openmm', 'integrator', 'steps_per_write']
        self.nblocks = (self.steps_per_tau // self.steps_per_write) + 1

        platform_name = config['west', 'openmm', 'platform', 'name'] or 'Reference'
        config_platform_properties = config['west', 'openmm', 'platform', 'properties'] or {}

        # Set up OpenMM
        with open(system_xml_file, 'r') as f:
            # NOTE: calling the system self.system causes a namespace collision in the propagator
            self.mmsystem = openmm.XmlSerializer.deserialize(f.read())

        with open(self.integrator_xml_file, 'r') as f:
            integrator = openmm.XmlSerializer.deserialize(f.read())

        self.platform = openmm.Platform.getPlatformByName(platform_name)
        self.platform_properties.update(config_platform_properties)

        self.temperature = integrator.getTemperature()

    @staticmethod
    def dist(x, y):
        return np.sqrt(np.sum((x-y)**2))

    @staticmethod
    def makepath(template, template_args=None,
                  expanduser=True, expandvars=True, abspath=False, realpath=False):
        template_args = template_args or {}
        path = template.format(**template_args)
        if expandvars: path = os.path.expandvars(path)
        if expanduser: path = os.path.expanduser(path)
        if realpath:   path = os.path.realpath(path)
        if abspath:    path = os.path.abspath(path)
        path = os.path.normpath(path)
        return path

    @staticmethod
    def mkdir_p(path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def get_pcoord(self, state):

        if isinstance(state, BasisState):
            coords = self.basis_coordinates.copy()
        elif isinstance(state, InitialState):
            template_args = {'initial_state': state}
            istate_data_ref = self.makepath(self.initial_state_ref_template, template_args)

            coords = np.loadtxt(istate_data_ref)
        else:
            raise TypeError('state must be BasisState or InitialState')

        state.pcoord = self.dist(coords[0,:], coords[1,:])

    def propagate(self, segments):

        platform_properties = self.platform_properties.copy()

        try:
            process_id = os.environ['WM_PROCESS_INDEX']
            platform_properties['OpenCLDeviceIndex'] = process_id
            platform_properties['CudaDeviceIndex'] = process_id
        except KeyError:
            pass

        with open(self.integrator_xml_file, 'r') as f:
            integrator = openmm.XmlSerializer.deserialize(f.read())
            integrator.setRandomNumberSeed(random.randint(0, 2**16))

        context = openmm.Context(self.mmsystem, integrator, self.platform, platform_properties)

        for segment in segments:
            starttime = time.time()

            # Set up arrays to hold trajectory data for pcoords, coordinates and velocities
            pcoords = np.empty((self.nblocks, 1))
            pcoords[0] = segment.pcoord[0]

            coordinates = np.empty((self.nblocks, self.mmsystem.getNumParticles(), 3))
            velocities = np.empty((self.nblocks, self.mmsystem.getNumParticles(), 3))

            # Get initial coordinates and velocities from restarts or initial state
            if segment.initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
                # Get restart data
                assert 'restart_coord' in segment.data
                assert 'restart_veloc' in segment.data

                coordinates[0] = segment.data['restart_coord']
                velocities[0] = segment.data['restart_veloc']

                initial_coords = units.Quantity(segment.data['restart_coord'], units.nanometer)
                initial_velocs = units.Quantity(segment.data['restart_veloc'], units.nanometer / units.picosecond)

                context.setPositions(initial_coords)
                context.setVelocities(initial_velocs)

                del segment.data['restart_coord']
                del segment.data['restart_veloc']

            elif segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ:
                initial_state = self.initial_states[segment.initial_state_id]

                assert initial_state.istate_type == InitialState.ISTATE_TYPE_GENERATED

                # Load coordinates coresponding to the initial state
                new_template_args = {'initial_state': initial_state}
                istate_data_ref = self.makepath(self.initial_state_ref_template, new_template_args)
                initial_coords = units.Quantity(np.loadtxt(istate_data_ref), units.angstrom)

                # Set up context for this segment
                context.setPositions(initial_coords)
                context.setVelocitiesToTemperature(self.temperature)

                state = context.getState(getPositions=True, getVelocities=True)
                coordinates[0] = state.getPositions(asNumpy=True)
                velocities[0] = state.getVelocities(asNumpy=True)

            # Run dynamics
            for istep in range(1, self.nblocks):
                integrator.step(self.steps_per_write)

                state = context.getState(getPositions=True, getVelocities=True)

                coordinates[istep] = state.getPositions(asNumpy=True)
                velocities[istep] = state.getVelocities(asNumpy=True)
                pcoords[istep] = 10.0 * self.dist(coordinates[istep,0,:], coordinates[istep,1,:])

            # Finalize segment trajectory
            segment.pcoord = pcoords[...].astype(pcoord_dtype)
            segment.data['coord'] = coordinates[...]
            segment.data['veloc'] = velocities[...]
            segment.status = Segment.SEG_STATUS_COMPLETE

            segment.walltime = time.time() - starttime

        return segments

    def gen_istate(self, basis_state, initial_state):
        '''Generate a new initial state from the given basis state.'''
        initial_coords = self.basis_coordinates.copy()
        initial_coords[0,0] = random.randrange(5, 16)

        new_template_args = {'initial_state': initial_state}
        istate_data_ref = self.makepath(self.initial_state_ref_template, new_template_args)
        self.mkdir_p(os.path.dirname(istate_data_ref))

        # Save coordinates of initial state as a text file
        # NOTE: this is ok for this example, but should be optimized for large systems
        np.savetxt(istate_data_ref, initial_coords)

        # Calculate pcoord for generated initial state
        pcoord = self.dist(initial_coords[0,:], initial_coords[1,:])
        initial_state.pcoord = np.array([pcoord], dtype=pcoord_dtype) 
        initial_state.istate_status = initial_state.ISTATE_STATUS_PREPARED

        return initial_state
