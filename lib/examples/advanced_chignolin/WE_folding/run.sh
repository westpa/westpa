#!/bin/bash

# Make sure environment is set
source env.sh

# Clean up
rm -f west.log

# Run w_run
$WEST_ROOT/bin/w_run --work-manager processes "$@" &> west.log
