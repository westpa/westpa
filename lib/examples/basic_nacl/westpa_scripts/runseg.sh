#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/bstate.pdb .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  ln -sv $WEST_PARENT_DATA_REF/seg.xml ./parent.xml
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  ln -sv $WEST_PARENT_DATA_REF ./parent.xml
fi

# Run the dynamics with OpenMM
python $WEST_SIM_ROOT/common_files/nacl_prod.py

#Calculate pcoord with MDAnalysis
python $WEST_SIM_ROOT/common_files/get_distance.py
cat dist.dat > $WEST_PCOORD_RETURN

# Clean up
rm -f parent.xml bstate.pdb dist.dat
