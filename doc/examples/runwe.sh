#!/bin/bash
#PBS -N we_test
#PBS -l walltime=4:00:00,mem=512mb,nodes=1:ppn=4

if [ -z "$TMPDIR" ] ; then
    TMPDIR=/scratch/mzwier/$$
    mkdir -p $TMPDIR
fi


SRCDIR=/scratch/mzwier/we_test
source $SRCDIR/env.sh
RUNSEG=$SRCDIR/runseg.sh
MULTIPLIER=8

cd $SIM_ROOT

time (while true; do
    echo
    date
    we_iter=$($WE cstatus -p we_iter) 
    if [ $we_iter -gt 2 ] ; then
        clean_iter=$[ $we_iter - 2 ]
        clean_dir=$SIM_ROOT/traj_segs/$clean_iter
        if ! [ -z "$clean_dir" ] ; then
            echo cleaning $clean_dir of restart data
            find $clean_dir -name \*.seg -o -name \*.trr | xargs rm -v
        else
            echo 'you are an idiot; check your script'
        fi
    elif [ $we_iter -gt 50 ] ; then 
        break
    fi
    $WE status
    while $WE cstatus -t segs_remaining; do
        for ((n=1; n<=$MULTIPLIER; n++)); do 
            for node in $(cat $PBS_NODEFILE); do
               ssh $node env SIM_ROOT=$SIM_ROOT TMPDIR=$TMPDIR bash $RUNSEG &
            done
        done
        wait
    done
    $WE we || exit 1
done)

