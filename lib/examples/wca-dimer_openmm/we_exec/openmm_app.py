import numpy as np

from simtk.openmm import app
from simtk.openmm import Vec3
import simtk.openmm.openmm as openmm
import simtk.unit as units

import wcadimer

from mdtraj.reporters import HDF5Reporter

import random
import argparse

from sys import stdout


def run(opts):
    system_xml_file = opts.system
    integrator_xml_file = opts.integrator
    coords_f = opts.coords
    platform_name = opts.platform
    deviceid = opts.device
    write_freq = opts.write_freq
    output = opts.output
    nsteps = opts.nsteps

    platform_properties = {'OpenCLPrecision': 'mixed',
                  'OpenCLPlatformIndex': '0',
                  'OpenCLDeviceIndex': '0',
                  'CudaPrecision': 'mixed',
                  'CudaDeviceIndex': '0',
                  'CpuThreads': '1'}

    platform_properties['CudaDeviceIndex'] = deviceid
    platform_properties['OpenCLDeviceIndex'] = deviceid
    
    with open(system_xml_file, 'r') as f:
        system = openmm.XmlSerializer.deserialize(f.read())

    with open(integrator_xml_file, 'r') as f:
        integrator = openmm.XmlSerializer.deserialize(f.read())
        integrator.setRandomNumberSeed(random.randint(0, 2**16))

    platform = openmm.Platform.getPlatformByName(platform_name)
    properties = {key: platform_properties[key] for key in platform_properties if key.lower().startswith(platform_name.lower())}
    if platform_name == 'CPU':
        properties = {'CpuThreads': '1'}

    print(properties)

    # Create dummy topology to satisfy Simulation object 
    topology = app.Topology()
    volume = wcadimer.natoms / wcadimer.density
    length = volume**(1.0/3.0)
    L = length.value_in_unit(units.nanometer)
    topology.setUnitCellDimensions(Vec3(L, L, L)*units.nanometer)

    simulation = app.Simulation(topology, system, 
                                integrator, platform, properties)

    init_data = np.load(coords_f)

    coords = units.Quantity(init_data['coord'], units.nanometer)
    simulation.context.setPositions(coords)

    if 'veloc' in init_data:
        velocs = units.Quantity(init_data['veloc'], units.nanometer / units.picosecond)
        simulation.context.setVelocities(velocs)
    else:
        simulation.context.setVelocitiesToTemperature(wcadimer.temperature)

    # Attach reporters
    simulation.reporters.append(HDF5Reporter(output + '.h5', write_freq, atomSubset=[0,1]))
    simulation.reporters.append(app.StateDataReporter(stdout, 20*write_freq, step=True, 
            potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
                speed=True, totalSteps=nsteps, separator='\t'))


    # Run segment
    simulation.step(nsteps)

    # Write restart data
    state = simulation.context.getState(getPositions=True, getVelocities=True)

    coords = state.getPositions(asNumpy=True)
    velocs = state.getVelocities(asNumpy=True)

    np.savez_compressed(output + '_restart.npz', coords, coord=coords, veloc=velocs)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coords', required=True, 
            help='Numpy npz file containing coordinates (and optionally velocities) to initialize simulation from')
    parser.add_argument('-s', '--system', required=True, help='OpenMM system xml file')
    parser.add_argument('-i', '--integrator', required=True, help='OpenMM integrator xml file')
    parser.add_argument('-p', '--platform', default='CUDA', 
            choices=['CUDA', 'OpenCL', 'Reference', 'CPU'], help='OpenMM platform name')
    parser.add_argument('-d', '--device', default='0', help='Device ID')
    parser.add_argument('-w', '--write_freq', required=True, type=int, help='write frequency for output files')
    parser.add_argument('-o', '--output', required=True, help='base name for output files')
    parser.add_argument('-n', '--nsteps', type=int, required=True, help='Number of steps to run')

    return parser

if __name__ == '__main__':
    parser = create_parser()
    opts = parser.parse_args()

    run(opts)
