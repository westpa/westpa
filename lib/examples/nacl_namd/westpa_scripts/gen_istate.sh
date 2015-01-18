#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)

# Generate random displacement from basis state
x=($(cat $WEST_BSTATE_DATA_REF | awk '{print $7;}'))
y=($(cat $WEST_BSTATE_DATA_REF | awk '{print $8;}'))
z=($(cat $WEST_BSTATE_DATA_REF | awk '{print $9;}'))
x[0]=$(echo "scale=4;${x[0]}-($WEST_RANDFLOAT/2)" | bc)
x[1]=$(echo "scale=4;${x[1]}+($WEST_RANDFLOAT/2)" | bc)

# Output to initial state
printf "ATOM      1  SOD SOD A   1     %7.3f %7.3f %7.3f\n" ${x[0]} ${y[0]} ${z[0]} \
  >> $WEST_ISTATE_DATA_REF
printf "ATOM      2  CLA CLA A   2     %7.3f %7.3f %7.3f\n" ${x[1]} ${y[1]} ${z[1]} \
  >> $WEST_ISTATE_DATA_REF
