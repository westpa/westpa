#!/bin/bash

cd $SIM_ROOT || exit 1

tar -r -f seg_logs.tar seg_logs/*.out && rm -f seg_logs/*.out


