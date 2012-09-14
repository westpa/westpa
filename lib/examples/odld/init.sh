#!/bin/bash

source env.sh

rm -f west.h5
$WEST_ROOT/bin/w_init --bstate initial,1.0 --tstate drift,10.01
