#!/usr/local/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
WCA fluid and WCA dimer systems.

DESCRIPTION

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

import numpy

import simtk
import simtk.unit as units
import simtk.openmm as openmm
    
#=============================================================================================
# CONSTANTS
#=============================================================================================

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

#=============================================================================================
# DEFAULT PARAMETERS
#=============================================================================================

natoms = 216                             # number of particles

# WCA fluid parameters (argon).
mass     = 39.9 * units.amu              # reference mass
sigma    = 3.4 * units.angstrom          # reference lengthscale
epsilon  = 120.0 * units.kelvin * kB     # reference energy

r_WCA    = 2.**(1./6.) * sigma           # interaction truncation range
tau = units.sqrt(sigma**2 * mass / epsilon) # characteristic timescale

# Simulation conditions.
temperature = 0.824 / (kB / epsilon)     # temperature
kT = kB * temperature                    # thermal energy
beta = 1.0 / kT                          # inverse temperature    
density  = 0.96 / sigma**3               # default density
stable_timestep = 0.001 * tau            # stable timestep
collision_rate = 1 / tau                 # collision rate for Langevin interator

# Dimer potential parameters.
h = 5.0 * kT                             # barrier height
r0 = r_WCA                               # compact state distance
w = 0.5 * r_WCA                          # extended state distance is r0 + 2*w

#=============================================================================================
# WCA Fluid
#=============================================================================================

def WCAFluid(N=natoms, density=density, mm=None, mass=mass, epsilon=epsilon, sigma=sigma):
    """
    Create a Weeks-Chandler-Andersen system.

    OPTIONAL ARGUMENTS

    N (int) - total number of atoms (default: 150)
    density (float) - N sigma^3 / V (default: 0.96)
    sigma
    epsilon

    """

    # Choose OpenMM package.
    if mm is None:
        mm = openmm

    # Create system
    system = mm.System()

    # Compute total system volume.
    volume = N / density
    
    # Make system cubic in dimension.
    length = volume**(1.0/3.0)
    # TODO: Can we change this to use tuples or 3x3 array?
    a = units.Quantity(numpy.array([1.0, 0.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    b = units.Quantity(numpy.array([0.0, 1.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    c = units.Quantity(numpy.array([0.0, 0.0, 1.0], numpy.float32), units.nanometer) * length/units.nanometer
    print("box edge length = %s" % str(length))
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Add particles to system.
    for n in range(N):
        system.addParticle(mass)
            
    # Create nonbonded force term implementing Kob-Andersen two-component Lennard-Jones interaction.
    energy_expression = '4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon'

    # Create force.
    force = mm.CustomNonbondedForce(energy_expression)

    # Set epsilon and sigma global parameters.
    force.addGlobalParameter('epsilon', epsilon)
    force.addGlobalParameter('sigma', sigma)

    # Add particles
    for n in range(N):
        force.addParticle([])    
    
    # Set periodic boundary conditions with cutoff.
    force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    print("setting cutoff distance to %s" % str(r_WCA))
    force.setCutoffDistance(r_WCA)    

    # Add nonbonded force term to the system.
    system.addForce(force)

    # Create initial coordinates using random positions.
    coordinates = units.Quantity(numpy.random.rand(N,3), units.nanometer) * (length / units.nanometer)
       
    # Return system and coordinates.
    return [system, coordinates]

#=============================================================================================
# WCA dimer plus fluid
#=============================================================================================

def WCADimer(N=natoms, density=density, mm=None, mass=mass, epsilon=epsilon, sigma=sigma, h=h, r0=r0, w=w):
    """
    Create a bistable bonded pair of particles (indices 0 and 1) optionally surrounded by a Weeks-Chandler-Andersen fluid.

    The bistable potential has form

    U(r) = h*(1-((r-r0-w)/w)^2)^2

    where r0 is the compact state separation, r0+2w is the extended state separation, and h is the barrier height.

    The WCA potential has form

    U(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ] + epsilon      (r < r*)
         = 0                                                       (r >= r*)

    where r* = 2^(1/6) sigma.

    OPTIONAL ARGUMENTS

    N (int) - total number of atoms (default: 2)
    density (float) - number density of particles (default: 0.96 / sigma**3)
    mass (simtk.unit.Quantity of mass) - particle mass (default: 39.948 amu)
    sigma (simtk.unit.Quantity of length) - Lennard-Jones sigma parameter (default: 0.3405 nm)
    epsilon (simtk.unit.Quantity of energy) - Lennard-Jones well depth (default: (119.8 Kelvin)*kB)
    h (simtk.unit.Quantity of energy) - bistable potential barrier height (default: ???)
    r0 (simtk.unit.Quantity of length) - bistable potential compact state separation (default: ???)
    w (simtk.unit.Quantity of length) - bistable potential extended state separation is r0+2*w (default: ???)

    """
    
    # Choose OpenMM package.
    if mm is None:
        mm = openmm

    # Compute cutoff for WCA fluid.
    r_WCA = 2.**(1./6.) * sigma # cutoff at minimum of potential

    # Create system
    system = mm.System()

    # Compute total system volume.
    volume = N / density
    
    # Make system cubic in dimension.
    length = volume**(1.0/3.0)
    a = units.Quantity(numpy.array([1.0, 0.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    b = units.Quantity(numpy.array([0.0, 1.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    c = units.Quantity(numpy.array([0.0, 0.0, 1.0], numpy.float32), units.nanometer) * length/units.nanometer
    print("box edge length = %s" % str(length))
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Add particles to system.
    for n in range(N):
        system.addParticle(mass)
            
    # WCA: Lennard-Jones truncated at minim and shifted so potential is zero at cutoff.
    energy_expression = '4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon'
    
    # Create force.
    force = mm.CustomNonbondedForce(energy_expression)

    # Set epsilon and sigma global parameters.
    force.addGlobalParameter('epsilon', epsilon)
    force.addGlobalParameter('sigma', sigma)

    # Add particles
    for n in range(N):
        force.addParticle([])

    # Add exclusion between bonded particles.
    force.addExclusion(0,1)
    
    # Set periodic boundary conditions with cutoff.
    if (N > 2):
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    else:
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    print("setting cutoff distance to %s" % str(r_WCA))
    force.setCutoffDistance(r_WCA)    

    # Add nonbonded force term to the system.
    system.addForce(force)

    # Add dimer potential to first two particles.
    dimer_force = openmm.CustomBondForce('h*(1-((r-r0-w)/w)^2)^2;')
    dimer_force.addGlobalParameter('h', h) # barrier height
    dimer_force.addGlobalParameter('r0', r0) # compact state separation
    dimer_force.addGlobalParameter('w', w) # second minimum is at r0 + 2*w
    dimer_force.addBond(0, 1, [])
    system.addForce(dimer_force)

    # Create initial coordinates using random positions.
    coordinates = units.Quantity(numpy.random.rand(N,3), units.nanometer) * (length / units.nanometer)
       
    # Reposition dimer particles at compact minimum.
    coordinates[0,:] *= 0.0
    coordinates[1,:] *= 0.0
    coordinates[1,0] = r0

    # Return system and coordinates.
    return [system, coordinates]

#=============================================================================================
# WCA dimer in vacuum
#=============================================================================================

def WCADimerVacuum(mm=None, mass=mass, epsilon=epsilon, sigma=sigma, h=h, r0=r0, w=w):
    """
    Create a bistable dimer.

    OPTIONAL ARGUMENTS

    """
    
    # Choose OpenMM package.
    if mm is None:
        mm = openmm

    # Create system
    system = mm.System()

    # Add particles to system.
    for n in range(2):
        system.addParticle(mass)

    # Add dimer potential to first two particles.
    dimer_force = openmm.CustomBondForce('h*(1-((r-r0-w)/w)^2)^2;')
    dimer_force.addGlobalParameter('h', h) # barrier height
    dimer_force.addGlobalParameter('r0', r0) # compact state separation
    dimer_force.addGlobalParameter('w', w) # second minimum is at r0 + 2*w
    dimer_force.addBond(0, 1, [])
    system.addForce(dimer_force)

    # Create initial coordinates using random positions.
    coordinates = units.Quantity(numpy.zeros([2,3], numpy.float64), units.nanometer)
       
    # Reposition dimer particles at compact minimum.
    coordinates[0,:] *= 0.0
    coordinates[1,:] *= 0.0
    coordinates[1,0] = r0

    # Return system and coordinates.
    return [system, coordinates]

