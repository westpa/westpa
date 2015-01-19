#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF
ln -sv $WEST_SIM_ROOT/gromacs_config/nacl.top .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/gromacs_config/md-continue.mdp > md.mdp
        ln -sv $WEST_PARENT_DATA_REF/seg.edr ./parent.edr
        ln -sv $WEST_PARENT_DATA_REF/seg.gro ./parent.gro
        ln -sv $WEST_PARENT_DATA_REF/seg.trr ./parent.trr
        $GROMPP -f md.mdp -c parent.gro -e parent.edr -p nacl.top \
          -t parent.trr -o seg.tpr -po md_out.mdp
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   appropriate basis or initial state
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/gromacs_config/md-genvel.mdp > md.mdp
        ln -sv $WEST_PARENT_DATA_REF ./parent.gro
        $GROMPP -f md.mdp -c parent.gro -p nacl.top \
          -o seg.tpr -po md_out.mdp
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
$MDRUN -s   seg.tpr -o seg.trr -c  seg.gro -e seg.edr \
       -cpo seg.cpt -g seg.log -nt 1

# Calculate progress coordinate
if [ ${G_DIST} ]; then
    # For GROMACS 4, use g_dist
    COMMAND="2 \n 3 \n"
    echo -e $COMMAND \
      | $G_DIST -f $WEST_CURRENT_SEG_DATA_REF/seg.trr -s seg.tpr -o dist.xvg
    cat dist.xvg | tail -n +23 | awk '{print $2*10;}' > $WEST_PCOORD_RETURN
elif [ ${GMX} ]; then
    # For GROMACS 5, use gmx distance
    $GMX distance -f $WEST_CURRENT_SEG_DATA_REF/seg.trr \
      -s $WEST_CURRENT_SEG_DATA_REF/seg.tpr -select 1 -oall dist.xvg
    cat dist.xvg | tail -n +16 | awk '{print $2*10;}' > $WEST_PCOORD_RETURN
fi

# Output coordinates
if [ ${WEST_COORD_RETURN} ]; then
    COMMAND="0 \n"
    if [ ${TRJCONV} ]; then
        # For GROMACS 4, use trjconv
        echo -e $COMMAND | $TRJCONV -f seg.trr -s seg.tpr -o seg.pdb
    elif [ ${GMX} ]; then
        # For GROMACS 5, use gmx trjconv
        echo -e $COMMAND | $GMX trjconv -f seg.trr -s seg.tpr -o seg.pdb
    fi
    cat $WEST_CURRENT_SEG_DATA_REF/seg.pdb | grep 'ATOM' \
      | awk '{print $6, $7, $8}' > $WEST_COORD_RETURN
fi

# Output log
if [ ${WEST_LOG_RETURN} ]; then
    cat $WEST_CURRENT_SEG_DATA_REF/seg.log \
      | awk '/Started mdrun/ {p=1}; p; /A V E R A G E S/ {p=0}' \
      > $WEST_LOG_RETURN
fi

# Clean up
rm -f dist.xvg md.mdp md_out.mdp nacl.top parent.gro parent.trr seg.cpt \
  seg.pdb seg.tpr
