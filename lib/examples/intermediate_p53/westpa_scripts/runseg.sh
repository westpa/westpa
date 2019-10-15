#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/P53.MDM2.prmtop .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi

$SANDER -O -i md.in   -p P53.MDM2.prmtop  -c parent.rst \
           -r seg.rst -x seg.nc      -o seg.log    -inf seg.nfo

COMMAND="         parm $WEST_SIM_ROOT/common_files/P53.MDM2.prmtop\n" 
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND reference $WEST_SIM_ROOT/bstates/P53.MDM2.rst\n"
COMMAND="$COMMAND rms p53-ca-rmsd :2-14@CA reference out ca-rmsd-p53.dat mass\n"
COMMAND="$COMMAND rms p53-heavy-rmsd :2-14&!@N,H,CA,HA,C,O reference out heavy-sc-rmsd-p53.dat mass\n"
COMMAND="$COMMAND multidihedral dihedralP53 phi psi resrange 2-14 out dihedral_p53.dat\n"
COMMAND="$COMMAND distance end-to-end :1 :15 out dist-end-to-end.dat\n"
COMMAND="$COMMAND go"

echo -e "$COMMAND" | cpptraj

paste <(cat ca-rmsd-p53.dat | tail -n +2 | awk {'print $2'}) <(cat dist-end-to-end.dat | tail -n +2 | awk {'print $2'})>$WEST_PCOORD_RETURN

cat heavy-sc-rmsd-p53.dat | tail -n +2 | awk {'print $2'} > $WEST_HEAVY_SC_RMSD_P53_RETURN

Everything is named DIHEDRAL_X
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $2'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $3'}) > $WEST_DIHEDRAL_2_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $4'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $5'}) > $WEST_DIHEDRAL_3_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $6'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $7'}) > $WEST_DIHEDRAL_4_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $8'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $9'}) > $WEST_DIHEDRAL_5_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $10'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $11'}) > $WEST_DIHEDRAL_6_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $12'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $13'}) > $WEST_DIHEDRAL_7_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $14'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $15'}) > $WEST_DIHEDRAL_8_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $16'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $17'}) > $WEST_DIHEDRAL_9_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $18'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $19'}) > $WEST_DIHEDRAL_10_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $20'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $21'}) > $WEST_DIHEDRAL_11_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $22'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $23'}) > $WEST_DIHEDRAL_12_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $24'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $25'}) > $WEST_DIHEDRAL_13_RETURN
paste <(cat dihedral_p53.dat | tail -n +2 | awk {'print $26'}) <(cat dihedral_p53.dat | tail -n +2 | awk {'print $27'}) > $WEST_DIHEDRAL_14_RETURN





# Clean up
rm -f $TEMP md.in parent.rst seg.nfo seg.pdb 

