#!/bin/bash
source env.sh
ps aux | grep w_run | grep -v grep

SFX=.d$$
mv seg_logs{,$SFX}
mv istates{,$SFX}
rm -Rf seg_logs$SFX istates$SFX & disown %1
rm -f system.h5 west.h5 seg_logs.tar
mkdir seg_logs istates

BSTATE_ARGS="--bstate unbound,1,nacl.inpcrd"
TSTATE_ARGS="--tstate bound,0.5"

$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 24 --work-manager=threads "$@"
