#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

DIST=_nacldist_$$.dat

function cleanup() {
    rm -f $DIST
}

trap cleanup EXIT

# Get progress coordinate
COMMAND="         trajin $WEST_STRUCT_DATA_REF\n"
COMMAND="$COMMAND distance na-cl :1@Cl- :2@Na+ out $DIST\n"
COMMAND="$COMMAND go\n"
echo -e $COMMAND | $CPPTRAJ  amber_config/nacl.prm
cat $DIST | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi
