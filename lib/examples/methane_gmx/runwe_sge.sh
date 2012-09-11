#!/bin/bash
#$ -N CH4bind
#$ -q compute48.q
#$ -l h_rt=36:00:00
#$ -pe 48way 48
#$ -V
#$ -cwd

# Set this to "1" if using multiple nodes; otherwise 0
MULTINODE=0

source /opt/modules/sgefix.sh
source env.sh

master=$(hostname)

# Switch depending on whether we are using multiple nodes
if [[ $MULTINODE == "1" ]] ; then
    # spawn workers
    for host in $(cat $TMPDIR/machines | uniq); do
        echo "starting WEMD node worker on $host"
        qrsh -inherit -V $host $PWD/node.sh --node -H $master </dev/null &> wemd-$JOB_ID-$host.log &
    done

    # spawn master
    echo "starting master on $master"
    $WEMD_ROOT/bin/w_run --master -H $master -n 0 &> wemd-$JOB_ID.log </dev/null &
    wait
else
    $WEMD_ROOT/bin/w_run --work-manager=threads &> wemd-$JOB_ID.log
fi

