#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEMD_SIM_ROOT

mkdir -pv $WEMD_CURRENT_SEG_DATA_REF || exit 1
cd $WEMD_CURRENT_SEG_DATA_REF || exit 1

TOP=methane.top
NDX=methane.ndx

# Set up the run
ln -sv $WEMD_SIM_ROOT/{$TOP,$NDX,*.mdp,*.itp} .

# if [ "$WEMD_PARENT_SEG_ID" -lt "0" ]; then
#     # this is a start/restart
#     ln -sv $WEMD_PARENT_SEG_DATA_REF/unbound.gro .
#     grompp -f md-genvel.mdp -c unbound.gro -p $TOP -o seg.tpr -n $NDX || exit 1
# else
#     ln -sv $WEMD_PARENT_SEG_DATA_REF/seg.gro ./parent.gro
#     ln -sv $WEMD_PARENT_SEG_DATA_REF/seg.trr ./parent.trr
#     grompp -f md-continue.mdp -c parent.gro -t parent.trr -p $TOP -o seg.tpr -n $NDX || exit 1
# fi

case $WEMD_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        ln -sv $WEMD_PARENT_DATA_REF/seg.gro ./parent.gro
        ln -sv $WEMD_PARENT_DATA_REF/seg.trr ./parent.trr
        grompp -f md-continue.mdp -c parent.gro -t parent.trr -p $TOP -o seg.tpr -n $NDX || exit 1
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEMD_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        ln -sv $WEMD_PARENT_DATA_REF ./initial.gro
        grompp -f md-genvel.mdp -c initial.gro -p $TOP -o seg.tpr -n $NDX || exit 1
    ;;

    *)
        echo "unknown init point type $WEMD_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac
    

# Propagate segment
mdrun -s seg.tpr -o seg.trr -x seg.xtc -c seg.gro \
      -e seg.edr -g seg.log -nt 1 \
      || exit 1

# Get progress coordinate
echo "7 8" | g_dist -f seg.xtc -s seg.tpr -n $NDX -xvg none || exit 1
awk '{print $2*10;}' < dist.xvg > $WEMD_PCOORD_RETURN || exit 1

# Get coordinates
if [ -n "$WEMD_COORD_RETURN" ] ; then
    echo 7 | g_traj -f seg.xtc -s seg.tpr -n $NDX -ox coord1.xvg -nopbc -fp -xvg none || exit 1
    echo 8 | g_traj -f seg.xtc -s seg.tpr -n $NDX -ox coord2.xvg -nopbc -fp -xvg none || exit 1
    paste coord1.xvg coord2.xvg > $WEMD_COORD_RETURN
fi

# Clean up
rm -f *.xvg *.itp *.mdp *.ndx *.top state.cpt

