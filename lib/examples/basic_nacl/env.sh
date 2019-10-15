#!/bin/bash

# Set up environment for westpa
export WEST_PYTHON=$(which python2.7)
source activate westpa-2017.10
export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)
