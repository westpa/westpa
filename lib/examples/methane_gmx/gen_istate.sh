#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEMD_SIM_ROOT

mkdir -p $(dirname $WEMD_ISTATE_DATA_REF)

# Simply copy in the basis state
#cp -v $WEMD_BSTATE_DATA_REF $WEMD_ISTATE_DATA_REF


./move_methane.py -o $WEMD_ISTATE_DATA_REF $(echo "$WEMD_RANDFLOAT*10.0" | bc) $WEMD_BSTATE_DATA_REF


