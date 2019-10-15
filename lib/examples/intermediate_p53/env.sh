#!/bin/bash

# Set up environment for dynamics
source /home/atb43/apps/amber18/amber.sh

# Set up environment for westpa
source /home/atb43/apps/westpa/westpa.sh
export WEST_PYTHON=$(which python2.7)
export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)

# Set runtime commands (this is said to be easier on the filesystem)
export SANDER=$(which sander)
export CPPTRAJ=$(which cpptraj)
