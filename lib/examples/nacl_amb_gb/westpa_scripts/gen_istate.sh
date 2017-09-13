#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)

# Generate random displacement from basis state
x=($(tail -n 1 $WEST_BSTATE_DATA_REF | awk '{print $1, $4;}'))
y=($(tail -n 1 $WEST_BSTATE_DATA_REF | awk '{print $2, $5;}'))
z=($(tail -n 1 $WEST_BSTATE_DATA_REF | awk '{print $3, $6;}'))
x[0]=$(echo "scale=7;${x[0]}-($WEST_RANDFLOAT/2)" | bc)
x[1]=$(echo "scale=7;${x[1]}+($WEST_RANDFLOAT/2)" | bc)

# Output to initial state
printf "sodium_chloride\n" >> $WEST_ISTATE_DATA_REF
printf "     2\n"          >> $WEST_ISTATE_DATA_REF
printf " %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n" ${x[0]} ${y[0]} ${z[0]} ${x[1]} ${y[1]} ${z[1]} \
  >> $WEST_ISTATE_DATA_REF
