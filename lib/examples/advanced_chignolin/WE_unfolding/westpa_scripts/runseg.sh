#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/chignolin.prmtop .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi

$SANDER -O -i md.in   -p chignolin.prmtop  -c parent.rst \
          -r seg.rst -x seg.nc      -o seg.log    -inf seg.nfo

TEMP=$(mktemp)
COMMAND="         parm chignolin.prmtop\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND reference $WEST_SIM_ROOT/common_files/chignolin.pdb\n"
COMMAND="$COMMAND rms ca-rmsd @CA reference out $TEMP mass\n"
COMMAND="$COMMAND go\n"

echo -e $COMMAND | $CPPTRAJ
cat $TEMP | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

cat $TEMP>pcoord.dat

# Clean up
rm -f $TEMP md.in parent.rst seg.nfo seg.pdb chignolin.prmtop
