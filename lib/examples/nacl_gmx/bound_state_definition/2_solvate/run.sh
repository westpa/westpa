gmx editconf \
  -f ../1_pdb2gmx/nacl_no_solvent_processed.gro\
  -c\
  -d 1.2\
  -bt octahedron\
  -o nacl_no_solvent_with_box.gro

gmx solvate \
  -cp nacl_no_solvent_with_box.gro\
  -cs spc216.gro\
  -p ../1_pdb2gmx/topol.top\
  -o nacl_solvated.gro
