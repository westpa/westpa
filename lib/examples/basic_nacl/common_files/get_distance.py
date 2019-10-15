import mdtraj
import numpy

traj = mdtraj.load('seg.dcd', top='bstate.pdb')
dist = mdtraj.compute_distances(traj, [[0,1]], periodic=True)
d_arr = numpy.asarray(dist)
d_arr = d_arr*10
numpy.savetxt("dist.dat", d_arr)
