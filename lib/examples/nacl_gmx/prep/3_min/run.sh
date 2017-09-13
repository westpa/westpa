#!/bin/sh
gmx grompp \
  -f minimize.mdp \
  -c ../2_solvate/nacl_solvated.gro \
  -p ../1_pdb2gmx/topol.top \
  -o 3_min.tpr

gmx mdrun \
  -deffnm 3_min\
  -mp ../1_pdb2gmx/topol.top\
  -nt 1\
  -v 

