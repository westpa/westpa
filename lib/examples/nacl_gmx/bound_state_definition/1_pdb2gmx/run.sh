#!/bin/sh

# Make the forcefield and symlink it here
cd ../../gromacs_config || exit 1
bash makeff.sh
cd ../bound_state_definition/1_pdb2gmx/ || exit 1
ln -s ../../gromacs_config/tip3p_ionsjc2008.ff .

gmx pdb2gmx \
  -f nacl_no_solvent.pdb \
  -ff tip3p_ionsjc2008 \
  -water tip3p \
  -o nacl_no_solvent_processed.gro
