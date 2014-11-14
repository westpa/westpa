#!/bin/bash

source env.sh

rm -f west.log
$WEST_ROOT/bin/w_run "$@" --n-workers=4 --work-manager=processes &>> west.log
