#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd        $WEST_CURRENT_SEG_DATA_REF || exit 1

PSF=nacl.psf
PRM=par_all27_prot_na.prm

# Set up the run
case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        INPUT_COOR=parent.coor
        INPUT_VEL=parent.vel
        CONF=md-continue.conf
        ln -sv $WEST_SIM_ROOT/namd_config/{$PSF,$PRM,$CONF} .
        ln -sv $WEST_PARENT_DATA_REF/seg.coor $INPUT_COOR
        ln -sv $WEST_PARENT_DATA_REF/seg.vel  $INPUT_VEL
        $NAMD $CONF > seg.log || exit 1
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        INPUT_COOR=initial.pdb
        CONF=md-genvel.conf
        ln -sv $WEST_SIM_ROOT/namd_config/{$PSF,$PRM,$CONF} .
        ln -sv $WEST_PARENT_DATA_REF $INPUT_COOR
        $NAMD $CONF > seg.log || exit 1
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

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
    if [ -n "$WEST_COORD_RETURN" ] ; then
        echo "${x[0]} ${y[0]} ${z[0]} ${x[1]} ${y[1]} ${z[1]}" >> $WEST_COORD_RETURN
    fi
done

# Clean up
rm -f $INPUT_COOR $INPUT_VEL *.conf *.dcd *.prm *.psf *.xsc [0-9]*
