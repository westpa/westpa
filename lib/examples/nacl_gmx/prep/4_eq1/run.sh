#!/bin/bash

ln -s ../1_pdb2gmx/posre.itp .

gmx grompp \
  -f 4_eq1.mdp \
  -c ../3_min/3_min.gro\
  -p ../1_pdb2gmx/topol.top\
  -o 4_eq1.tpr

gmx mdrun \
  -deffnm 4_eq1\
  -mp ../1_pdb2gmx/topol.top\
  -v
