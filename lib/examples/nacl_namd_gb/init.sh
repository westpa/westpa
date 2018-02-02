#!/bin/bash
source env.sh
ps aux | grep w_run | grep -v grep

SFX=.d$$
mv      traj_segs{,$SFX}
mv      seg_logs{,$SFX}
mv      istates{,$SFX}
rm -Rf  traj_segs$SFX seg_logs$SFX istates$SFX & disown %1
rm -f   system.h5 west.h5 seg_logs.tar
mkdir   seg_logs traj_segs istates

BSTATE_ARGS="--bstate-file bstates/bstates.txt"
TSTATE_ARGS="--tstate bound,1.8 --tstate unbound,16.9"

$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 24 \
  --work-manager=threads "$@"
