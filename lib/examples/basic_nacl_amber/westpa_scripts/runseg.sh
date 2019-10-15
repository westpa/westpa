#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/nacl.parm7 .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi

$PMEMD -O -i md.in   -p nacl.parm7  -c parent.rst \
          -r seg.rst -x seg.nc      -o seg.log    -inf seg.nfo

TEMP=$(mktemp)
COMMAND="         parm nacl.parm7\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND autoimage fixed Na+ \n"
COMMAND="$COMMAND distance na-cl :1@Na+ :2@Cl- out $TEMP\n"
COMMAND="$COMMAND go\n"

echo -e $COMMAND | $CPPTRAJ
cat $TEMP | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

if [ ${WEST_COORD_RETURN} ]; then
  COMMAND="         parm nacl.parm7\n"
  COMMAND="$COMMAND trajin  $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
  COMMAND="$COMMAND strip :WAT \n"
  COMMAND="$COMMAND autoimage fixed Na+ \n"
  COMMAND="$COMMAND trajout $WEST_CURRENT_SEG_DATA_REF/seg.pdb\n"
  COMMAND="$COMMAND go\n"
  echo -e $COMMAND | $CPPTRAJ 
  cat $WEST_CURRENT_SEG_DATA_REF/seg.pdb | grep 'ATOM' \
    | awk '{print $6, $7, $8}' > $WEST_COORD_RETURN
fi

# Clean up
rm -f $TEMP md.in parent.rst seg.nfo seg.pdb nacl.parm7
