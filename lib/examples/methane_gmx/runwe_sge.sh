#!/bin/bash
#$ -N CH4bind
#$ -V
#$ -cwd

source env.sh || exit 1

cd $WEST_SIM_ROOT
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$JOB_ID.json

# start server
$WEST_ROOT/bin/w_run --work-manager=zmq --n-workers=0 --zmq-mode=server --zmq-info=$SERVER_INFO &> west-$JOB_ID.log &

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
while read -r machine ncores _j _k ; do 
   qrsh -inherit -V $machine $PWD/node.sh --work-manager=zmq --zmq-mode=client --n-workers=$ncores --zmq-info=$SERVER_INFO &> west-$JOB_ID-$machine.log &
done < $PE_HOSTFILE

wait
