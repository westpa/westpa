#!/usr/bin/env python
import MDAnalysis
import MDAnalysis.analysis.distances

universe = MDAnalysis.Universe('../1_psf/nacl.psf', '2_min.dcd')
nacl = universe.select_atoms('resname SOD or resname CLA')
outfile = open('dist.dat','w+')
for ts in universe.trajectory:
    dist = MDAnalysis.analysis.distances.self_distance_array(nacl.atoms.positions)[0]
    outfile.write('{:.03f}\n'.format(dist))

print("Final distance: {:03f} Angstroms".format(dist))
