#!/usr/bin/env python
import MDAnalysis

# Load the trajectories for the current segment and the parent segment.
# We need all timepoints from the current segment, and the final timepoint
# from the parent segment (which is where the current segment starts).
parent_universe = MDAnalysis.Universe('structure.psf', 'parent.dcd')
current_universe = MDAnalysis.Universe('structure.psf', 'seg.dcd')


with MDAnalysis.Writer('seg.pdb', multiframe=True, bonds=None, n_atoms=2) as PDB:
    # Go to the last timepoint of the parent trajectory and output a pdb.
    parent_nacl = parent_universe.select_atoms('resname SOD or resname CLA')
    parent_universe.trajectory[-1]
    PDB.write(parent_nacl.atoms)

    # Now write out pdb frames for each timepoint of the current segment.
    nacl = current_universe.select_atoms('resname SOD or resname CLA')
    for ts in current_universe.trajectory:
        PDB.write(nacl.atoms)

