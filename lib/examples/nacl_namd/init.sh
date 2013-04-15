#!/bin/bash
source env.sh
ps aux | grep w_run | grep -v grep
pkill -9 -f w_run

SFX=.d$$
mv      traj_segs{,$SFX}
mv      seg_logs{,$SFX}
mv      istates{,$SFX}
rm -Rf  traj_segs$SFX seg_logs$SFX istates$SFX & disown %1
rm -f   system.h5 west.h5 seg_logs.tar
mkdir   seg_logs traj_segs istates

BSTATE_ARGS="--bstate-file bstates/bstates.txt"
TSTATE_ARGS="--tstate bound,2.8"

$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 48 --work-manager=threads "$@"
