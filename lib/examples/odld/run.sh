#!/bin/bash

source env.sh

rm -f wemd.log
$WEMD_ROOT/bin/w_run --verbose &> wemd.log
