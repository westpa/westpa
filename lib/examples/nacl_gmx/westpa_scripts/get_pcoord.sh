#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

DIST=_nacldist_$$.xvg

function cleanup() {
    rm -f $DIST
}

trap cleanup EXIT

# Get progress coordinate
if [ ${G_DIST} ]; then
    # For GROMACS 4, use g_dist
    COMMAND="2 \n 3 \n"
    echo -e $COMMAND | $G_DIST -f $WEST_STRUCT_DATA_REF \
      -s $WEST_SIM_ROOT/gromacs_config/nacl.tpr -o $DIST
    cat $DIST | tail -n +23 | awk '{print $2*10;}' > $WEST_PCOORD_RETURN
elif [ ${GMX} ]; then
    # For GROMACS 5, use gmx distance
    $GMX distance -f $WEST_STRUCT_DATA_REF  -s $WEST_STRUCT_DATA_REF \
      -oall $DIST -select 1
    cat $DIST | tail -n +16 | awk '{print $2*10;}' > $WEST_PCOORD_RETURN
fi

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi
