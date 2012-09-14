#!/bin/bash
#$ -N CH4bind
#$ -q compute48.q
#$ -l h_rt=36:00:00
#$ -pe 8way 16
#$ -V
#$ -cwd

source /opt/modules/sgefix.sh
source env.sh || exit 1

cd $WEST_SIM_ROOT
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$JOB_ID.json

# start server
$WEST_ROOT/bin/w_run --wm-work-manager=zmq --wm-n-workers=0 --wm-zmq-mode=server --wm-zmq-server-info=$SERVER_INFO &> west-$JOB_ID.log &

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
   qrsh -inherit -V $machine $PWD/node-ltc1.sh --wm-work-manager=zmq --wm-zmq-mode=client --wm-n-workers=$ncores --wm-zmq-server-info=$SERVER_INFO &> west-$JOB_ID-$machine.log &
done < $PE_HOSTFILE

wait
