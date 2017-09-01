#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)
ln -s $WEST_BSTATE_DATA_REF $WEST_ISTATE_DATA_REF

