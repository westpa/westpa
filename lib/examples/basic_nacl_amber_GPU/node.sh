#!/bin/bash -l

set -x

#cd $WEST_SIM_ROOT
cd $1; shift
source env.sh
export WEST_JOBID=$1; shift
export SLURM_NODENAME=$1; shift
# export LOCAL=/local/$WEST_JOBID
echo "starting WEST client processes on: "; hostname
echo "current directory is $PWD"
echo "environment is: "
env | sort

$WEST_ROOT/bin/w_run "$@" &> west-$SLURM_NODENAME-node.log

echo "Shutting down.  Hopefully this was on purpose?"
