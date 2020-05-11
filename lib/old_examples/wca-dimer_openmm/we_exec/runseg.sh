#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1

SYS=system.xml
INTEGRATOR=integrator.xml
OMM_APP=openmm_app.py
GET_PCOORD=calc_pcoord.py

# Set up the run
ln -sv $WEST_SIM_ROOT/{$SYS,$INTEGRATOR,$OMM_APP,$GET_PCOORD} .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        ln -sv $WEST_PARENT_DATA_REF/seg_restart.npz ./initial.npz
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        ln -sv $WEST_PARENT_DATA_REF ./initial.npz
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac
    

# Propagate segment
$WEST_PYTHON $OMM_APP -c initial.npz -s $SYS \
       -i $INTEGRATOR -d $WM_PROCESS_INDEX \
       -p CUDA -w 250 -o seg -n 250 || exit 1

# Get progress coordinate
$WEST_PYTHON $GET_PCOORD -f initial.npz -o d1.dat || exit 1
$WEST_PYTHON $GET_PCOORD -f seg.h5 -o d2.dat || exit 1
cat d1.dat d2.dat > $WEST_PCOORD_RETURN
cp $WEST_PCOORD_RETURN pcoord.dat


# Clean up
rm -f *.xml *.py initial.npz

