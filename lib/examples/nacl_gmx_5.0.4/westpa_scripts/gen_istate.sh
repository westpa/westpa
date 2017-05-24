#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)

# Generate random displacement from basis state
x=($(cat $WEST_BSTATE_DATA_REF | grep '[NA,CL]' | awk '{print $4}'))
y=($(cat $WEST_BSTATE_DATA_REF | grep '[NA,CL]' | awk '{print $5}'))
z=($(cat $WEST_BSTATE_DATA_REF | grep '[NA,CL]' | awk '{print $6}'))
x[0]=$(echo "scale=3;((${x[0]}*10)-($WEST_RANDFLOAT/2))/10" | bc)
x[1]=$(echo "scale=3;((${x[1]}*10)+($WEST_RANDFLOAT/2))/10" | bc)

# Output to initial state
printf "sodium_chloride\n" >  $WEST_ISTATE_DATA_REF
printf "    2\n"           >> $WEST_ISTATE_DATA_REF
printf "    1Na      NA    1%8.3f%8.3f%8.3f\n" ${x[0]} ${y[0]} ${z[0]} \
  >> $WEST_ISTATE_DATA_REF
printf "    2Cl      CL    2%8.3f%8.3f%8.3f\n" ${x[1]} ${y[1]} ${z[1]} \
  >> $WEST_ISTATE_DATA_REF
printf "   0.99000   0.00000   0.00000\n" >> $WEST_ISTATE_DATA_REF
