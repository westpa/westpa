#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

OUT=dist_$$.dat

function cleanup() {
    rm -f $OUT
}

trap cleanup EXIT

# Get progress coordinate
$WEST_PYTHON calc_pcoord.py -f $WEST_STRUCT_DATA_REF -o $WEST_PCOORD_RETURN || exit 1
#cat < $OUT > $WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi


