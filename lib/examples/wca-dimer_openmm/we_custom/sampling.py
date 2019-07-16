import simtk
import simtk.unit as units
import simtk.openmm as openmm


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
