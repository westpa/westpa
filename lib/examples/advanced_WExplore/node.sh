#!/bin/bash -l

set -x

cd $SLURM_SUBMIT_DIR
(source env.sh
cd $WEST_SIM_ROOT

$WEST_ROOT/bin/w_run "$@" ) &> west-node.log
