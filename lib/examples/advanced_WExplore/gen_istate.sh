#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)

# Simply copy in the basis state
cp -v $WEST_BSTATE_DATA_REF.gro $WEST_ISTATE_DATA_REF.gro
cp -v $WEST_BSTATE_DATA_REF.trr $WEST_ISTATE_DATA_REF.trr
cp -v $WEST_BSTATE_DATA_REF.edr $WEST_ISTATE_DATA_REF.edr

