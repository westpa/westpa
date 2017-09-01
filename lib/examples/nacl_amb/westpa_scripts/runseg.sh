#!/bin/bash
#
# runseg.sh
#
# WESTPA runs this script for each trajectory segment. WESTPA supplies
# environment variables that are unique to each segment, such as:
#
#   WEST_CURRENT_SEG_DATA_REF: A path to where the current trajectory segment's
#       data will be stored. This will become "WEST_PARENT_DATA_REF" for any
#       child segments that spawn from this segment
#   WEST_PARENT_DATA_REF: A path to a file or directory containing data for the
#       parent segment.
#   WEST_CURRENT_SEG_INITPOINT_TYPE: Specifies whether this segment is starting
#       anew, or if this segment continues from where another segment left off.
#   WEST_RAND16: A random integer
#
# This script has the following three jobs:
#  1. Create a directory for the current trajectory segment, and set up the
#     directory for running pmemd/sander 
#  2. Run the dynamics
#  3. Calculate the progress coordinates and return data to WESTPA


# If we are running in debug mode, then output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

######################## Set up for running the dynamics #######################

# Set up the directory where data for this segment will be stored.
cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

# Make a symbolic link to the topology file. This is not unique to each segment.
ln -sv $WEST_SIM_ROOT/amber_config/nacl.parm7 .

# Either continue an existing tractory, or start a new trajectory. In the 
# latter case, we need to do a couple things differently, such as generating
# velocities.
#
# First, take care of the case that this segment is a continuation of another
# segment.  WESTPA provides the environment variable 
# $WEST_CURRENT_SEG_INITPOINT_TYPE, and we check its value.
if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  # The weighted ensemble algorithm requires that dynamics are stochastic.
  # We'll use the "sed" command to replace the string "RAND" with a randomly
  # generated seed.
  sed "s/RAND/$WEST_RAND16/g" \
      $WEST_SIM_ROOT/amber_config/md.in > md.in

  # This trajectory segment will start off where its parent segment left off.
  # The "ln" command makes symbolic links to the parent segment's rst file. 
  # This is preferable to copying the files, since it doesn't
  # require writing all the data again.
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst

# Now take care of the case that the trajectory is starting anew.
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  # Again, we'll use the "sed" command to replace the string "RAND" with a 
  # randomly generated seed.
  sed "s/RAND/$WEST_RAND16/g" \
      $WEST_SIM_ROOT/amber_config/md.in > md.in
  # For a new segment, we only need to make a symbolic link to the .rst file.
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi

############################## Run the dynamics ################################
# Propagate segment using pmemd (or sander)
$PMEMD -O -i md.in   -p nacl.parm7  -c parent.rst \
          -r seg.rst -x seg.nc      -o seg.log    -inf seg.nfo

########################## Calculate and return data ###########################

# Calculate progress coordinate
# Note that it is very important to include the parent.rst file!
# WESTPA assumes that the first timepoint of a given iteration is the same as
# the final timepoint of the previous iteration.
TEMP=$(mktemp)
COMMAND="         parm nacl.parm7\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND distance na-cl :1@Na+ :2@Cl- out $TEMP\n"
COMMAND="$COMMAND go\n"
echo -e $COMMAND | $CPPTRAJ
cat $TEMP | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

# Output coordinates
if [ ${WEST_COORD_RETURN} ]; then
  COMMAND="         parm nacl.parm7\n"
  COMMAND="$COMMAND trajin  $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
  COMMAND="$COMMAND trajin  $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
  COMMAND="$COMMAND strip :WAT \n"
  COMMAND="$COMMAND autoimage fixed Na+ \n"
  COMMAND="$COMMAND trajout $WEST_CURRENT_SEG_DATA_REF/seg.pdb\n"
  COMMAND="$COMMAND go\n"
  echo -e $COMMAND | $CPPTRAJ 
  cat $WEST_CURRENT_SEG_DATA_REF/seg.pdb | grep 'ATOM' \
    | awk '{print $6, $7, $8}' > $WEST_COORD_RETURN
fi

# Clean up
rm -f $TEMP md.in parent.rst seg.nfo seg.pdb nacl.parm7
