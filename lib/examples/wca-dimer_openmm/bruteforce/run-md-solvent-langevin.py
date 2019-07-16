

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Simulation of WCA dimer in dense WCA solvent using Langevin Dynamics

DESCRIPTION

Modified from the original code developed by John Choder 
The original unmodified code is available at:
    https://simtk.org/home/ncmc

Modifications made by Joshua L. Adelman


COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

TODO

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import errno
import sys
import math
import copy
import time

import numpy
import numpy as np

import simtk
import simtk.unit as units
import simtk.openmm as openmm
    
#import Scientific.IO.NetCDF as netcdf # for netcdf interface in Scientific
import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

import wcadimer

#=============================================================================================
# SUBROUTINES
#=============================================================================================
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def norm(n01):
    return n01.unit * numpy.sqrt(numpy.dot(n01/n01.unit, n01/n01.unit))


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


def equilibrate_langevin(system, timestep, collision_rate, temperature, sqrt_kT_over_m, coordinates, platform):
    nsteps = 5000

    print("Equilibrating for %.3f ps..." % ((nsteps * timestep) / units.picoseconds))
    
    # Create integrator and context.
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)

    # Set coordinates.
    context.setPositions(coordinates)

    # Set Maxwell-Boltzmann velocities
    velocities = sqrt_kT_over_m * numpy.random.standard_normal(size=sqrt_kT_over_m.shape)
    context.setVelocities(velocities)

    # Equilibrate.
    integrator.step(nsteps)

    # Compute energy
    print("Computing energy.")
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()
    print("potential energy: %.3f kcal/mol" % (potential_energy / units.kilocalories_per_mole))

    # Get coordinates.
    state = context.getState(getPositions=True, getVelocities=True)    
    coordinates = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)    
    box_vectors = state.getPeriodicBoxVectors()
    system.setDefaultPeriodicBoxVectors(*box_vectors)    

    print("Computing energy again.")
    context.setPositions(coordinates)
    context.setVelocities(velocities)        
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()
    print("potential energy: %.3f kcal/mol" % (potential_energy / units.kilocalories_per_mole))
    
    total_energy = compute_energy(context, coordinates, velocities)    

    return [coordinates, velocities]


def compute_energy(context, positions, velocities):
    """
    Compute total energy for positions and velocities.
    """
    # Set positions and velocities.
    context.setPositions(positions)
    context.setVelocities(velocities)
    # Compute total energy.
    state = context.getState(getEnergy=True)
    total_energy = state.getPotentialEnergy() + state.getKineticEnergy()

    #print "potential energy: %.3f kcal/mol" % (state.getPotentialEnergy() / units.kilocalories_per_mole)
    #print "kinetic   energy: %.3f kcal/mol" % (state.getKineticEnergy() / units.kilocalories_per_mole)    
    
    return total_energy

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    mkdir_p('data')

    # PARAMETERS
    netcdf_filename = 'data/md-solvent-langevin.nc'


    verbose = False
    
    # WCA fluid parameters (argon).
    mass     = wcadimer.mass
    sigma    = wcadimer.sigma
    epsilon  = wcadimer.epsilon
    r_WCA    = wcadimer.r_WCA
    r0       = wcadimer.r0
    h        = wcadimer.h
    w        = wcadimer.w
    
    # Compute characteristic timescale.
    tau = wcadimer.tau
    print("tau = %.3f ps" % (tau / units.picoseconds))

    # Compute timestep.
    equilibrate_timestep = 2 * wcadimer.stable_timestep
    timestep = equilibrate_timestep
    print("equilibrate timestep = %.1f fs, switch timestep = %.1f fs" % (equilibrate_timestep / units.femtoseconds, timestep / units.femtoseconds))
    print("timestep: %f" % (timestep / units.femtoseconds))

    # Set temperature, pressure, and collision rate for stochastic thermostats.
    temperature = wcadimer.temperature
    print("temperature = %.1f K" % (temperature / units.kelvin))
    kT = wcadimer.kT
    beta = 1.0 / kT # inverse temperature    
    collision_rate = 1.0 / tau # collision rate for Langevin integrator
    print('collision_rate: ', collision_rate)

    niterations = 1000000  # number of work samples to collect
    nsteps = 250 # number of steps per interation
    deviceid = 0
    print('nsteps: ', nsteps)

    # Create system.     
    [system, coordinates] = wcadimer.WCADimer()

    # Form vectors of masses and sqrt(kT/m) for force propagation and velocity randomization.
    print("Creating masses array...")
    nparticles = system.getNumParticles()
    masses = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.amu)
    for particle_index in range(nparticles):
        masses[particle_index,:] = system.getParticleMass(particle_index)
    sqrt_kT_over_m = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.nanometers / units.picosecond)
    for particle_index in range(nparticles):
        sqrt_kT_over_m[particle_index,:] = units.sqrt(kT / masses[particle_index,0]) # standard deviation of velocity distribution for each coordinate for this atom

    # List all available platforms
    print("Available platforms:")
    for platform_index in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(platform_index)
        print("%5d %s" % (platform_index, platform.getName()))
    print("")

    # Select platform.
    #platform = openmm.Platform.getPlatformByName("CPU")
    platform = openmm.Platform.getPlatformByName("CUDA")
    min_platform = openmm.Platform.getPlatformByName("Reference")

    for prop in platform.getPropertyNames():
        print(prop, platform.getPropertyDefaultValue(prop))

    platform.setPropertyDefaultValue('CudaDeviceIndex', '%d' % deviceid)         
    platform.setPropertyDefaultValue('CudaPrecision', 'mixed')
    #platform.setPropertyDefaultValue('CpuThreads', '1')

    # Initialize netcdf file.
    if not os.path.exists(netcdf_filename):
        # Open NetCDF file for writing
        ncfile = netcdf.Dataset(netcdf_filename, 'w') # for netCDF4
        
        # Create dimensions.
        ncfile.createDimension('iteration', 0) # unlimited number of iterations

        # Create variables.
        ncfile.createVariable('distance', 'd', ('iteration',))
        
        # Force sync to disk to avoid data loss.
        ncfile.sync()

        # Minimize.
        print("Minimizing energy...")
        coordinates = minimize(min_platform, system, coordinates)
    
        # Equilibrate.
        print("Equilibrating...")
        coordinates, velocities = equilibrate_langevin(system, equilibrate_timestep, collision_rate, temperature, sqrt_kT_over_m, coordinates, platform)
        
        # Write initial configuration.
        np.save('data/restart.npy', coordinates[:,:] / units.angstroms)

        ncfile.variables['distance'][0] = norm(coordinates[1,:] - coordinates[0,:]) / units.angstroms
        ncfile.sync()        
        iteration = 1
    else:
        # Open NetCDF file for reading.
        ncfile = netcdf.Dataset(netcdf_filename, 'a') # for netCDF4

        # Read iteration and coordinates.
        iteration = ncfile.variables['distance'][:].size
        coordinates = units.Quantity(np.load('data/restart.npy'), units.angstroms)

    # Continue
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    context = openmm.Context(system, integrator, platform)
    context.setPositions(coordinates)
    context.setVelocitiesToTemperature(temperature)

    state = context.getState(getPositions=True)

    while (iteration < niterations):
        #print 'coordinates: ', coordinates[0,:] / units.angstrom
        initial_time = time.time()
        
        # Generate a new configuration.
        initial_distance = norm(coordinates[1,:] - coordinates[0,:])

        integrator.step(nsteps)

        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        
        final_distance = norm(coordinates[1,:] - coordinates[0,:])            
        
        ncfile.variables['distance'][iteration] = final_distance / units.angstroms
        ncfile.sync()

        np.save('data/restart.npy', coordinates[:,:] / units.angstroms)

        # Debug.
        final_time = time.time()
        elapsed_time = final_time - initial_time
        if iteration % 100 == 0:
            print("iteration %5d / %5d" % (iteration, niterations))
            print("Dynamics %.1f A -> %.1f A (barrier at %.1f A)" % (initial_distance / units.angstroms, final_distance / units.angstroms, (r0+w)/units.angstroms))
            print("%12.3f s elapsed" % elapsed_time)

        sys.stdout.flush()
        sys.stderr.flush()

        # Increment iteration counter.
        iteration += 1

    # Close netcdf file.
    ncfile.close()
