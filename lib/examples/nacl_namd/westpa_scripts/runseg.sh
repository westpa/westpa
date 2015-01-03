#!/bin/bash

set -ex

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1
ln -sv $WEST_SIM_ROOT/namd_config/par_all27_prot_na.prm .
ln -sv $WEST_SIM_ROOT/namd_config/nacl.psf              .
ln -sv $WEST_SIM_ROOT/namd_config/nacl.pdb              .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/namd_config/md-continue.conf \
          > md.conf
        ln -sv $WEST_PARENT_DATA_REF/seg.coor              ./parent.coor
        ln -sv $WEST_PARENT_DATA_REF/seg.vel               ./parent.vel
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        # A  appropriate basis or initial state
        sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/namd_config/md-genvel.conf \
          > md.conf
        ln -sv $WEST_PARENT_DATA_REF                     ./parent.coor
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
$NAMD md.conf > seg.log || exit 1

# Calculate progress coordinate and output coordinates
frames=($(ls [0-9][0-9][0-9][0-9].coor))
for frame in ${frames[@]}; do
    x=($(cat $frame | head -n 3 | tail -n 2 | cut -c 31-38))
    y=($(cat $frame | head -n 3 | tail -n 2 | cut -c 39-46))
    z=($(cat $frame | head -n 3 | tail -n 2 | cut -c 47-54))
    dx=$(echo "(${x[0]})-(${x[1]})" | bc)
    dy=$(echo "(${y[0]})-(${y[1]})" | bc)
    dz=$(echo "(${z[0]})-(${z[1]})" | bc)
    echo "sqrt (($dx*$dx)+($dy*$dy)+($dz*$dz))" | bc >> $WEST_PCOORD_RETURN
    [ ${WEST_COORD_RETURN} ] &&
        echo "${x[0]} ${y[0]} ${z[0]} ${x[1]} ${y[1]} ${z[1]}" >> $WEST_COORD_RETURN
done

# Output log
[ ${WEST_LOG_RETURN} ] &&
    cat seg.log | grep ^ENERGY | tail -n +2 | cut -c 8- > $WEST_LOG_RETURN

# Clean up
rm -f [0-9]* md.conf nacl.pdb nacl.psf par_all27_prot_na.prm parent.coor
