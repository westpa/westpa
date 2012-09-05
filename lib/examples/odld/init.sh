#!/bin/bash

source env.sh

rm -f wemd.h5
$WEMD_ROOT/bin/w_init --bstate initial,1.0 --tstate drift,10.01
