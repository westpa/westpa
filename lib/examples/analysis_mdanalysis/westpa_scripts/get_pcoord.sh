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
cd $WEST_SIM_ROOT/common_files

# Set the arguments for rmsd.py and call the script to calculate initial progress
# coordinate.
# Arguments:
#   ref: path to initial state coordinate file.
#   top: path to topology file.
#   mob: path to trajectory file.
#   for: we are evaluating a basis/initial state, so for = 'RESTRT'
$WEST_SIM_ROOT/rmsd.py \
    --ref $WEST_SIM_ROOT/bstates/P53.MDM2.rst \
    --top $WEST_SIM_ROOT/amber_config/P53.MDM2.prmtop \
    --mob $WEST_SIM_ROOT/bstates/P53.MDM2.rst \
    --for RESTRT \

cat rmsd.dat > $WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
