#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT || exit 1

ITER=$(printf "%06d" $WEST_CURRENT_ITER)
tar -cf seg_logs/$ITER.tar seg_logs/$ITER-*.log
rm  -f  seg_logs/$ITER-*.log

BACKUPINTERVAL=100
CHECK=$(python -c "print($WEST_CURRENT_ITER%$BACKUPINTERVAL)")
if [ "${CHECK}" == "0" ]; then
  cp ${WEST_SIM_ROOT}/west.h5 ${WEST_SIM_ROOT}/west.h5.backup1
  mv ${WEST_SIM_ROOT}/west.h5.backup1 ${WEST_SIM_ROOT}/west.h5.backup2
fi

