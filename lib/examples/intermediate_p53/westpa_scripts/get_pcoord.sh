#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT/common_files

COMMAND="           parm $WEST_SIM_ROOT/common_files/P53.MDM2.prmtop\n" 
COMMAND="$COMMAND trajin $WEST_STRUCT_DATA_REF\n"
COMMAND="$COMMAND reference $WEST_SIM_ROOT/bstates/P53.MDM2.rst\n"
COMMAND="$COMMAND rms p53-ca-rmsd :2-14@CA reference out ca-rmsd-p53.dat mass\n"
COMMAND="$COMMAND rms p53-heavy-rmsd :2-14&!@N,H,CA,HA,C,O reference out heavy-sc-rmsd-p53.dat mass\n"
COMMAND="$COMMAND multidihedral dihedralP53 phi psi resrange 2-14 out dihedral_p53.dat\n"
COMMAND="$COMMAND distance end-to-end :1 :15 out dist-end-to-end.dat\n"
COMMAND="$COMMAND go"

echo -e "$COMMAND" | cpptraj

paste <(cat ca-rmsd-p53.dat | tail -n +2 | awk {'print $2'}) <(cat dist-end-to-end.dat | tail -n +2 | awk {'print $2'})>$WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi



