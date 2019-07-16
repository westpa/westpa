#!/bin/sh
#
# env.sh
#
# This script defines environment variables that are used by other shell
# scripts, both when setting up the simulation and when running the simulation.
#

################################## NAMD ########################################

# Set an environment variables for namd2. This is good practice for running
# WESTPA simulations on clusters, as repeatedly calling "namd2" rather than the 
# absolute path to the binaries (e.g., "/usr/local/namd/namd2") can take a long 
# time and be harder on the filesystem.  For this tutorial, this step is not 
# necessary but nonetheless demonstrates good practice.
export NAMD=$(which namd2)

############################## Python and WESTPA ###############################
# Next inform WESTPA what python it should use.  
export WEST_PYTHON=$(which python3)

# Check to make sure that the environment variable WEST_ROOT is set. 
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.
if [ -z "$WEST_ROOT" ]; then
  echo "The environment variable WEST_ROOT is not set."
  echo "Try running 'source westpa.sh' from the WESTPA installation directory"
  exit 1
fi

# Explicitly name our simulation root directory.  Similar to the statement 
# above, we check if the variable is not set.  If the variable is not set,
# the we set it 
if [ -z "$WEST_SIM_ROOT" ]; then
  export WEST_SIM_ROOT="$PWD"
fi

# Set the simulation name.  Whereas "WEST_SIM_ROOT" gives us the entire 
# absolute path to the simulation directory, running the "basename" command
# will give us only the last part of that path (the directory name).
export SIM_NAME=$(basename $WEST_SIM_ROOT)
