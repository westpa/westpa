#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        ln -sv $WEST_PARENT_DATA_REF/odld.log ./parent.log
        ln -sv $WEST_PARENT_DATA_REF/odld.rst ./odld.crd
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   appropriate basis or initial state
        ln -sv $WEST_PARENT_DATA_REF ./odld.crd
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
python $WEST_SIM_ROOT/odld.py >& odld.log

# Calculate progress coordinate
cat pcoords.dat > $WEST_PCOORD_RETURN

cat displacements.dat > $WEST_DISPLACEMENT_RETURN

