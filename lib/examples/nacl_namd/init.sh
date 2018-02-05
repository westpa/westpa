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

# Copy over the equilibrated conformation from ./prep to bstates/unbound,
# including coordinates, velocities, and box information
mkdir bstates/unbound
cp prep/4_eq2/4_eq2.coor bstates/unbound/seg.coor
cp prep/4_eq2/4_eq2.dcd  bstates/unbound/seg.dcd
cp prep/4_eq2/4_eq2.vel  bstates/unbound/seg.vel
cp prep/4_eq2/4_eq2.xsc  bstates/unbound/seg.xsc

cp prep/1_psf/nacl.psf   namd_config/nacl.psf
cp prep/1_psf/nacl.pdb   namd_config/nacl.pdb

# Make sure the CHARMM force field is in the right spot:
if [ ! -d namd_config/toppar ]; then
  echo "CHARMM force field not found! Download the CHARMM36 force field, untar "
  echo "it, and place the 'toppar' directory in namd_config."
  exit 1
fi

# Process the CHARMM force field file to be compatibile with NAMD
cat namd_config/toppar/toppar_water_ions.str \
  | grep -v -e "^set" \
  | grep -v -e "^if" \
  | grep -v -e "^WRNLEV"\
  | grep -v -e "^BOMLEV"\
  | grep -v -e "^SOD\s\s\s\sO"\
  > namd_config/toppar/toppar_water_ions_for_namd.str

# Define the arguments for the basis states (used for generating initial 
# states; in this case we only have one), and target states (used for
# knowing when to recycle trajectories). In this example, we recycle
# trajectories as they reach the bound state; we focus on sampling  
# the binding process (and not the unbinding process).

BSTATE_ARGS="--bstate-file bstates/bstates.txt"
TSTATE_ARGS="--tstate bound,1.0"

# Initialize the simulation, creating the main WESTPA data file (west.h5)
# The "$@" lets us take any arguments that were passed to init.sh at the
# command line and pass them along to w_init.
$WEST_ROOT/bin/w_init \
  $BSTATE_ARGS $TSTATE_ARGS \
  --segs-per-state 5 \
  --work-manager=threads "$@"
