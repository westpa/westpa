#!/bin/bash
#SBATCH --job-name=GPU
#SBATCH --output=slurm.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cluster=invest
#SBATCH --partition=lchong
#SBATCH --gres=gpu:2
#SBATCH --time=24:00:00

set -x
cd $SLURM_SUBMIT_DIR
source env.sh || exit 1

env | sort

cd $WEST_SIM_ROOT
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json

# start server
$WEST_ROOT/bin/w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west-$SLURM_JOBID.log &

# wait on host info file up to one minute
for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

# start clients, with the proper number of cores on each

for node in $(scontrol show hostname $SLURM_NODELIST); do
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node --work-manager=zmq --zmq-mode=client --zmq-read-host-info=$SERVER_INFO --zmq-comm-mode=tcp &
done


wait
