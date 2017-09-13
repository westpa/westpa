#!/usr/bin/env python
import MDAnalysis
import MDAnalysis.analysis.distances

# Load the trajectories for the current segment and the parent segment.
# We need all timepoints from the current segment, and the final timepoint
# from the parent segment (which is where the current segment starts)
parent_universe = MDAnalysis.Universe('structure.psf', 'parent.dcd')
current_universe = MDAnalysis.Universe('structure.psf', 'seg.dcd')

# Go to the last timepoint of the parent trajectory and calculate the distance
# between the ions
parent_nacl = parent_universe.select_atoms('resname SOD or resname CLA')
parent_universe.trajectory[-1]
dist = MDAnalysis.analysis.distances.self_distance_array(parent_nacl.atoms.positions, 
                                                         box=parent_universe.dimensions)[0]
print("{:.03f}".format(dist))

# Now calculate the distance between the ions for each timepoint of the current
# trajectory.
nacl = current_universe.select_atoms('resname SOD or resname CLA')
for ts in current_universe.trajectory:
    dist = MDAnalysis.analysis.distances.self_distance_array(nacl.atoms.positions,
                                                             box=current_universe.dimensions)[0]
    print("{:.03f}".format(dist))
