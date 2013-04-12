#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

PTRAJ=_ptraj_$$.tmp
DIST=_nacldist_$$.dat

function cleanup() {
    rm -f $PTRAJ $DIST
}

trap cleanup EXIT

# Get progress coordinate
echo "trajin $WEST_STRUCT_DATA_REF" >> $PTRAJ
echo "distance na-cl :1@Cl- :2@Na+ out $DIST" >> $PTRAJ
ptraj nacl.prmtop $PTRAJ
awk '{print $2}' < $DIST > $WEST_PCOORD_RETURN || exit 1

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi


