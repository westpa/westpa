from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys

pdb = PDBFile('bstate.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.loadState('parent.xml')
simulation.reporters.append(StateDataReporter('seg.log', 100, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True)) 
simulation.saveState('seg.xml')
simulation.reporters.append(DCDReporter('seg.dcd', 20)) 
simulation.step(1000)
