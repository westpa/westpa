#!/bin/bash

cd $WEST_SIM_ROOT || exit 1

rm -f seg_logs/*.log

# save mdcrd and rst files from each segment; toss everything else
EXTS="edr gro trr xtc"
if [[ $WEST_CURRENT_ITER -gt 2 ]] ; then
    del_iter=$[ $WEST_CURRENT_ITER - 2]
    del_iter_dir=traj_segs/$(printf "%06d" $del_iter)
    pushd $del_iter_dir &> /dev/null || exit 1
#        for seg_id in *; do
#            pushd $seg_id &> /dev/null || exit 1
#                # protect files to save
#                for ext in $EXTS; do
#                    mv seg.$ext .seg.$ext
#                done
#                # delete everything else
#                rm *
#                # deprotect files to save
#                for ext in $EXTS; do
#                    mv .seg.$ext seg.$ext
#                done
#            popd &> /dev/null || exit 1
#        done
    popd &> /dev/null || exit 1
    #echo -n "tar traj_segs: "; time tar -c $del_iter_dir > $del_iter_dir.tar && rm -Rf $del_iter_dir
    rm -Rf $del_iter_dir
fi
