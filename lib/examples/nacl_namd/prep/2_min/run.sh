#!/bin/bash

# Make a pdb that specifies the harmonic restraints on Na+ and Cl-;
# This is the same as the original pdb file, with the Bfactor column
# edited to reflect the force constants for each restraint
cat ../1_psf/nacl.pdb \
  | sed -e "/SOD\|CLA/ s/\s0.00\s/ 5.00 /" \
  > restraints.pdb

namd2 2_min.conf > 2_min.log

