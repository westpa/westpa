#!/bin/sh
gmx grompp \
  -f minimize.mdp \
  -c ../2_solvate/nacl_solvated.gro \
  -p ../1_pdb2gmx/topol.top \
  -o minimization.tpr

gmx mdrun \
  -deffnm minimization\
  -mp ../1_pdb2gmx/topol.top\
  -v \
  -nt 1
