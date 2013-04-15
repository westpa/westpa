#!/bin/bash
if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi
cd $WEST_SIM_ROOT

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1


# Set up the run
ln -sv $WEST_SIM_ROOT/{nacl.prmtop,ptraj.continue} .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        ln -sv $WEST_SIM_ROOT/md-continue.in ./md.in
        ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        cp $WEST_SIM_ROOT/md-genvel.in ./md.in
        cp $WEST_PARENT_DATA_REF  ./parent.rst
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
     sander -O -i md.in -p nacl.prmtop -c parent.rst -r seg.rst -x seg.mdcrd -o seg.out || exit 1


# Get progress coordinate
ptraj nacl.prmtop ptraj.continue
awk '{print $2}' < nacldist.dat > $WEST_PCOORD_RETURN

#Get coordinates
tail < seg.mdcrd > $WEST_COORD_RETURN

# Clean up
rm -f nacl* md.in ptraj.continue

