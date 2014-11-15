import simtk.openmm.openmm as openmm
import wcadimer


system, _ = wcadimer.WCADimer()

temperature = wcadimer.temperature
collision_rate = wcadimer.collision_rate
timestep = 2.0 * wcadimer.stable_timestep
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

# Serialize openmm objects
with open('system.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(system))

with open('integrator.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(integrator))
