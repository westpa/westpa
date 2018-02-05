#!/bin/bash
#
# init.sh
#
# Initialize the WESTPA simulation, creating initial states (istates) and the
# main WESTPA data file, west.h5. 
#
# If you run this script after starting the simulation, the data you generated
# will be erased!
#

source env.sh

# Make sure that seg_logs (log files for each westpa segment), traj_segs (data
# from each trajectory segment), and istates (initial states for starting new
# trajectories) directories exist and are empty. 
rm -rf traj_segs seg_logs istates west.h5 
mkdir   seg_logs traj_segs istates

# Define the arguments for the basis states (used for generating initial 
# states; in this case we only have one), and target states (used for
# knowing when to recycle trajectories).
BSTATE_ARGS="--bstate-file bstates/bstates.txt"
TSTATE_ARGS="--tstate bound,1.8 --tstate unbound,16.9"

# Initialize the simulation, creating the main WESTPA data file (west.h5)
# The "$@" lets us take any arguments that were passed to init.sh at the
# command line and pass them along to w_init.
$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 24 \
  --work-manager=threads "$@"
