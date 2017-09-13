#!/bin/bash

ln -s ../1_pdb2gmx/posre.itp .

gmx grompp \
  -f 5_eq2.mdp \
  -c ../4_eq1/4_eq1.gro\
  -t ../4_eq1/4_eq1.cpt \
  -p ../1_pdb2gmx/topol.top\
  -o 5_eq2.tpr

gmx mdrun \
  -deffnm 5_eq2\
  -mp ../1_pdb2gmx/topol.top\
  -v
