import simtk.openmm.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as units


amber_prmtop = app.AmberPrmtopFile('nacl.prmtop')
amber_inpcrd = app.AmberInpcrdFile('nacl.inpcrd')

system = amber_prmtop.createSystem(implicitSolvent=app.HCT,)
integrator = openmm.LangevinIntegrator(300*units.kelvin, 0.5/units.picoseconds, 2.0*units.femtoseconds)

# Serialize openmm objects
with open('system.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(system))

with open('integrator.xml', 'w') as f:
    f.write(openmm.XmlSerializer.serialize(integrator))
