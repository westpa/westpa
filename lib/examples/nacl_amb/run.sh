#!/bin/bash

source env.sh

rm -f west.log
$WEST_ROOT/bin/w_run "$@" &> west.log
