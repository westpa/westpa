#!/bin/bash

source env.sh

rm -f west.h5
BSTATES="--bstate initial,1.0"
#TSTATES="--tstate drift,10.01 --tstate bound,1.3"
#TSTATES="--tstate bound,1.2"
$WEST_ROOT/bin/w_init $BSTATES $TSTATES "$@"
