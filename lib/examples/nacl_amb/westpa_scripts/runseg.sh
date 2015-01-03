#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1
ln -sv $WEST_SIM_ROOT/amber_config/nacl.prm .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/amber_config/md-continue.in \
          > md.in
        ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   appropriate basis or initial state
        sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/amber_config/md-genvel.in \
          > md.in
        ln -sv $WEST_PARENT_DATA_REF ./parent.rst
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
$PMEMD -O -i md.in   -p nacl.prm  -c parent.rst \
          -r seg.rst -x seg.crd   -o seg.log    -inf seg.nfo || exit 1


# Calculate progress coordinate
COMMAND="         trajin $WEST_CURRENT_SEG_DATA_REF/seg.crd\n"
COMMAND="$COMMAND distance na-cl :1@Cl- :2@Na+ out dist.dat\n"
COMMAND="$COMMAND go\n"
echo -e $COMMAND | $CPPTRAJ nacl.prm
cat dist.dat | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

# Output coordinates and log
[ ${WEST_COORD_RETURN} ] &&
    echo $WEST_CURRENT_SEG_DATA_REF/seg.crd > $WEST_COORD_RETURN
[ ${WEST_LOG_RETURN} ] &&
    echo $WEST_CURRENT_SEG_DATA_REF/seg.log > $WEST_LOG_RETURN

# Clean up
rm -f dist.dat md.in nacl.prm parent.rst seg.nfo
