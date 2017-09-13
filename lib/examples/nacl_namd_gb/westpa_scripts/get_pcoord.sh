#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Get progress coordinate
x=($(cat $WEST_STRUCT_DATA_REF | cut -c 31-38))
y=($(cat $WEST_STRUCT_DATA_REF | cut -c 39-46))
z=($(cat $WEST_STRUCT_DATA_REF | cut -c 47-54))
dx=$(echo "(${x[0]})-(${x[1]})" | bc)
dy=$(echo "(${y[0]})-(${y[1]})" | bc)
dz=$(echo "(${z[0]})-(${z[1]})" | bc)
echo "sqrt (($dx*$dx)+($dy*$dy)+($dz*$dz))" | bc > $WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi
