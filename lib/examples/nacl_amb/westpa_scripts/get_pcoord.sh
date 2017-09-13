#!/bin/bash
#
# get_pcoord.sh
#
# This script is run when calculating initial progress coordinates for new 
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT

# Make a temporary file in which to store output from cpptraj
DIST=$(mktemp)


# Load the restart (.rst) file into cpptraj and calculate the distance between
# the Na+ and Cl- ions. Here, $WEST_STUCT_DATA_REF indicates a particular 
# $WEST_ISTATE_DATA_REF, as defined by gen_istate.sh
COMMAND="           parm $WEST_SIM_ROOT/amber_config/nacl.parm7 \n"
COMMAND="${COMMAND} trajin $WEST_STRUCT_DATA_REF \n"
COMMAND="${COMMAND} distance na-cl :1@Na+ :2@Cl- out $DIST \n"
COMMAND="${COMMAND} go"
echo -e "${COMMAND}" | $CPPTRAJ

# Pipe the relevant part of the output file (the distance) to $WEST_PCOORD_RETURN
cat $DIST | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

# Remove the temporary file to clean up
rm $DIST

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
