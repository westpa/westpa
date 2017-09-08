#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF
ln -sv $WEST_SIM_ROOT/namd_config/par_all27_prot_na.prm .
ln -sv $WEST_SIM_ROOT/namd_config/nacl.psf              .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/namd_config/md-continue.conf > md.conf
        ln -sv $WEST_PARENT_DATA_REF/seg.coor ./parent.coor
        ln -sv $WEST_PARENT_DATA_REF/seg.vel  ./parent.vel
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   appropriate basis or initial state
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/namd_config/md-genvel.conf > md.conf
        ln -sv $WEST_PARENT_DATA_REF ./parent.coor
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
$NAMD md.conf > seg.log

# Calculate progress coordinate and output coordinates
frames=(parent.coor $(ls [0-9][0-9][0-9][0-9].coor))
for frame in ${frames[@]}; do
    x=($(cat $frame | head -n 3 | tail -n 2 | cut -c 31-38))
    y=($(cat $frame | head -n 3 | tail -n 2 | cut -c 39-46))
    z=($(cat $frame | head -n 3 | tail -n 2 | cut -c 47-54))
    dx=$(echo "(${x[0]})-(${x[1]})" | bc)
    dy=$(echo "(${y[0]})-(${y[1]})" | bc)
    dz=$(echo "(${z[0]})-(${z[1]})" | bc)
    echo "sqrt (($dx*$dx)+($dy*$dy)+($dz*$dz))" | bc >> $WEST_PCOORD_RETURN
    if [ ${WEST_COORD_RETURN} ]; then
        echo "${x[0]} ${y[0]} ${z[0]} ${x[1]} ${y[1]} ${z[1]}" \
          >> $WEST_COORD_RETURN
    fi
done

# Output log
if [ ${WEST_LOG_RETURN} ]; then
    cat seg.log | grep ^ENERGY | cut -c 8- > $WEST_LOG_RETURN
fi

# Clean up
rm -f [0-9]* md.conf nacl.psf par_all27_prot_na.prm parent.coor
