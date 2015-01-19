#!/bin/bash

cd $WEST_SIM_ROOT
source env.sh

set -x
$WEST_ROOT/bin/w_run "$@"
