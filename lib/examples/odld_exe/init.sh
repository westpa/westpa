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

BSTATES="--bstate initial,1.0"
$WEST_ROOT/bin/w_init $BSTATES $TSTATES "$@"
