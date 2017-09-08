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

# Make a temporary file to hold the output from "gmx distance"
DIST=$(mktemp)

# Calculate the distance between Na+ and Cl-
$GMX distance -f $WEST_STRUCT_DATA_REF -s $WEST_STRUCT_DATA_REF -oall $DIST \
  -select 1

# Return the distance between Na+ and Cl- WESTPA.
#                   |keep the last line only                
#                   |                |-Take the second column
#                   |                |  |-multiply the distance by 10 to
#                   V                V  V convert nm to Angstroms.
cat $WEST_STRUCT_DATA_REF > $WEST_SIM_ROOT/test.txt
cat $DIST.xvg >> $WEST_SIM_ROOT/test.txt
cat $DIST.xvg | tail -n 1 | awk '{print $2*10;}' > $WEST_PCOORD_RETURN

# Remove the temporary file to clean up.
rm $DIST.xvg

