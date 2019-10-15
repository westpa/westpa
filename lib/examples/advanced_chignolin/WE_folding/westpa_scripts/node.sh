#!/bin/bash
#
# node.sh

cd $WEST_SIM_ROOT
source env.sh

set -x
$WEST_ROOT/bin/w_run "$@"
