#!/bin/sh

# Check for gmx
#GMX=$(which gmx)
#if [ ! -n "${GMX}"]; then
#  echo "Gromacs not found! Make sure your Gromacs installation is in your path."
#  exit 1
#fi
#GROMACS_DIR=$(basename ${GMX} /bin/gmx)
#cp ${GROMACS_DIR}/share/gromacs/top/amber03.ff/tip3p.itp tip3p_jc2008.ff/
ln -s ../../gromacs_config/tip3p_ionsjc2008.ff tip3p_ionsjc2008.ff

gmx pdb2gmx \
  -f nacl_no_solvent.pdb \
  -ff tip3p_ionsjc2008 \
  -water tip3p \
  -o nacl_no_solvent_processed.gro
