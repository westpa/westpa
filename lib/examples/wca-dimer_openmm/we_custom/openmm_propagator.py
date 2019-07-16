import os
import errno
import random
import time
import glob
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

pcoord_dtype = np.float32


class OpenMMPropagator(WESTPropagator):
    def __init__(self, rc=None):
        super(OpenMMPropagator, self).__init__(rc)

        self.pcoord_dtype = pcoord_dtype
        self.pcoord_ndim = 1

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
        self.basis_state_ref_template = config['west','data','data_refs','basis_state']

        system_xml_file = config['west', 'openmm', 'system', 'file']
        self.integrator_xml_file = config['west', 'openmm', 'integrator', 'file']

        self.steps_per_tau = config['west', 'openmm', 'integrator', 'steps_per_tau']
        self.steps_per_write = config['west', 'openmm', 'integrator', 'steps_per_write']
        self.nblocks = (self.steps_per_tau // self.steps_per_write) + 1

        platform_name = config['west', 'openmm', 'platform', 'name'] or 'Reference'

        # Set up OpenMM
        with open(system_xml_file, 'r') as f:
            # NOTE: calling the system self.system causes a namespace collision in the propagator
            self.mmsystem = openmm.XmlSerializer.deserialize(f.read())

        with open(self.integrator_xml_file, 'r') as f:
            integrator = openmm.XmlSerializer.deserialize(f.read())

        self.platform = openmm.Platform.getPlatformByName(platform_name)

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
            template_args = {'basis_state': state}
            bstate_data_ref = self.makepath(self.basis_state_ref_template, template_args)
            print(bstate_data_ref)

            coords = 10.0 * np.load(bstate_data_ref)
        elif isinstance(state, InitialState):
            template_args = {'initial_state': state}
            istate_data_ref = self.makepath(self.initial_state_ref_template, template_args)

            coords = 10.0 * np.load(istate_data_ref)
        else:
            raise TypeError('state must be BasisState or InitialState')

        state.pcoord = self.dist(coords[0,:], coords[1,:])
        print(state.pcoord)

    @staticmethod
    def load_parent_data(n_iter):
        restart_files = glob.glob('traj_segs/iter_{:06d}_*.npz'.format(n_iter))
        parent_coords = {}
        parent_velocs = {}
        for rf in restart_files:
            data = np.load(rf)
            for si, seg_id in enumerate(data['seg_ids']):
                parent_coords[seg_id] = data['coords'][si]
                parent_velocs[seg_id] = data['velocs'][si]

        return parent_coords, parent_velocs

    def propagate(self, segments):

        platform_properties = {key: value for key, value in self.platform_properties.items() if key.startswith(self.platform.getName())}

        try:
            process_id = os.environ['WM_PROCESS_INDEX']
            if self.platform.getName() == 'OpenCL':
                platform_properties['OpenCLDeviceIndex'] = process_id
            elif self.platform.getName() == 'CUDA':
                platform_properties['CudaDeviceIndex'] = process_id
            elif self.platform.getName() == 'CPU':
                platform_properties['CpuThreads'] = '1'
        except KeyError:
            process_id = 0

        with open(self.integrator_xml_file, 'r') as f:
            integrator = openmm.XmlSerializer.deserialize(f.read())
            integrator.setRandomNumberSeed(random.randint(0, 2**16))

        context = openmm.Context(self.mmsystem, integrator, self.platform, platform_properties)

        if segments[0].n_iter > 1:
            parent_coords, parent_velocs = self.load_parent_data(segments[0].n_iter - 1)

        block_coordinates = np.empty((len(segments), self.mmsystem.getNumParticles(), 3))
        block_velocities = np.empty((len(segments), self.mmsystem.getNumParticles(), 3))
        block_seg_ids = np.empty(len(segments), dtype=np.int)

        for si, segment in enumerate(segments):
            starttime = time.time()

            # Set up arrays to hold trajectory data for pcoords, coordinates and velocities
            pcoords = np.empty((self.nblocks, 1))
            pcoords[0] = segment.pcoord[0]

            coordinates = np.empty((self.nblocks, self.mmsystem.getNumParticles(), 3))
            velocities = np.empty((self.nblocks, self.mmsystem.getNumParticles(), 3))

            # Get initial coordinates and velocities from restarts or initial state
            if segment.initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
                # Get restart data
                coordinates[0] = parent_coords[segment.parent_id]
                velocities[0] = parent_velocs[segment.parent_id]

                initial_coords = units.Quantity(parent_coords[segment.parent_id], units.nanometer)
                initial_velocs = units.Quantity(parent_velocs[segment.parent_id], units.nanometer / units.picosecond)

                context.setPositions(initial_coords)
                context.setVelocities(initial_velocs)

            elif segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ:
                initial_state_id = segment.initial_state_id
                basis_state_id = self.initial_states[initial_state_id].basis_state_id

                assert basis_state_id in [0,1]
                if basis_state_id == 0:
                    tag = '_a'
                else:
                    tag = '_b'

                basis_fname = os.path.join(os.environ['WEST_SIM_ROOT'], 'bstates', 'init_coords{}.npy'.format(tag))
                initial_coords = units.Quantity(np.load(basis_fname), units.nanometer)

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

                # Check for system blowing up
                assert pcoords[istep] < 12.0, 'pcoord dist: {}'.format(pcoords[istep])
                assert coordinates[istep].max() < 20.0, 'max coord: {}'.format(coordinates[istep].max())

            # Finalize segment trajectory
            segment.pcoord = pcoords[...].astype(pcoord_dtype)
            segment.status = Segment.SEG_STATUS_COMPLETE

            block_coordinates[si] = coordinates[-1]
            block_velocities[si] = velocities[-1]
            block_seg_ids[si] = segment.seg_id

            segment.walltime = time.time() - starttime

        np.savez_compressed('traj_segs/iter_{:06d}_{:06d}.npz'.format(segments[0].n_iter, block_seg_ids[0]), 
                coords=block_coordinates,
                velocs=block_velocities,
                seg_ids=block_seg_ids)

        return segments

