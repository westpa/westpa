import os
import numpy as np
import argparse

import simtk
import simtk.unit as units
import simtk.openmm as openmm

import wcadimer

def minimize(platform, system, positions):
    # Create a Context.
    timestep = 1.0 * units.femtoseconds
    integrator = openmm.VerletIntegrator(timestep)
    context = openmm.Context(system, integrator, platform)
    # Set coordinates.
    context.setPositions(positions)
    # Compute initial energy.
    state = context.getState(getEnergy=True)
    initial_potential = state.getPotentialEnergy()
    print("initial potential: %12.3f kcal/mol" % (initial_potential / units.kilocalories_per_mole))
    # Minimize.
    openmm.LocalEnergyMinimizer.minimize(context)
    # Compute final energy.
    state = context.getState(getEnergy=True, getPositions=True)
    final_potential = state.getPotentialEnergy()
    positions = state.getPositions(asNumpy=True)
    # Report
    print("final potential  : %12.3f kcal/mol" % (final_potential / units.kilocalories_per_mole))

    return positions


def norm(n01):
    return n01.unit * np.sqrt(np.dot(n01/n01.unit, n01/n01.unit))


def run(platform_name, deviceid, two):
    not_obs = [True, True]

    system, coordinates = wcadimer.WCADimer()
    print("Time step: ", (wcadimer.stable_timestep * 2.0).in_units_of(units.femtoseconds))

    # Minimization
    platform = openmm.Platform.getPlatformByName('Reference')

    print('Minimizing energy...')
    coordinates = minimize(platform, system, coordinates)
    print('Separation distance: {}'.format(norm(coordinates[1,:] - coordinates[0,:]) / units.angstroms))

    print('Equilibrating...')
    platform = openmm.Platform.getPlatformByName(platform_name)

    if platform_name == 'CUDA':
        platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % deviceid)
        platform.setPropertyDefaultValue('CudaPrecision', 'mixed')
    elif platform_name == 'OpenCL':
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % deviceid)
        platform.setPropertyDefaultValue('OpenCLPrecision', 'mixed')

    integrator = openmm.LangevinIntegrator(wcadimer.temperature, 
            wcadimer.collision_rate, 
            2.0 * wcadimer.stable_timestep)

    context = openmm.Context(system, integrator, platform)
    context.setPositions(coordinates)
    context.setVelocitiesToTemperature(wcadimer.temperature)

    if two:
        while not_obs[0] or not_obs[1]:
            integrator.step(5000)

            state = context.getState(getPositions=True)
            coordinates = state.getPositions(asNumpy=True)
            sep_dist = norm(coordinates[1,:] - coordinates[0,:]) / units.angstroms
            print('Separation distance: {}'.format(sep_dist))
            if sep_dist < 5.7:
                not_obs[0] = False
                tag = '_a'
                sep_dist_a = sep_dist
            else:
                not_obs[1] = False
                tag = '_b'
                sep_dist_b = sep_dist

            if not os.path.isdir('bstates'):
                os.makedirs('bstates')

            np.save('bstates/init_coords{}.npy'.format(tag), coordinates / units.nanometers)

        print(sep_dist_a, sep_dist_b)

    else:
        integrator.step(5000)

        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        print('Separation distance: {}'.format(norm(coordinates[1,:] - coordinates[0,:]) / units.angstroms))

        if not os.path.isdir('bstates'):
            os.makedirs('bstates')
            np.save('bstates/init_coords.npy', coordinates / units.nanometers)


if __name__ == '__main__':
    available_platforms = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]

    print('Available platforms:')
    for platform_index, platform in enumerate(available_platforms):
        print("%5d %s" % (platform_index, platform))
    print()

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--platform', choices=available_platforms,
            default='CUDA',
            help='Platform to use for equilibration')
    parser.add_argument('-d', '--deviceid', type=int, default=0, help='Device id for equilibration if using a GPU')
    parser.add_argument('--two', action="store_true", default=False, help='Create initial states in both basins')

    args = parser.parse_args()

    run(args.platform, args.deviceid, args.two)
