#!/bin/bash

[ ${WEST_SIM_ROOT} ] ||
    exit 1
[ ! -d $WEST_SIM_ROOT/traj_segs ] &&
    exit 1

cd $WEST_SIM_ROOT/traj_segs

ITERS=($(ls | grep '^[0-9][0-9][0-9][0-9][0-9][0-9]$'))
ITERS=("${ITERS[@]:0:${#ITERS[@]}-1}")
for ITER in ${ITERS[@]}; do
    [ ! -f $ITER.tar ] &&
        tar -cvf $ITER.tar $ITER
done

