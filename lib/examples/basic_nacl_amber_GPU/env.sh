# Modify to taste

# Load Modules
module purge
module load intel/2017.3.196 amber/18

module unload python
module load cuda/8.0.44
#module load intel-mpi
export CUDA_HOME=/opt/packages/cuda/8.0
export LD_LIBRARY_PATH=/opt/packages/cuda/8.0/lib64:$LD_LIBRARY_PATH

# This is our local scratch, where we'll store files during the dynamics.
export NODELOC=$LOCAL
export USE_LOCAL_SCRATCH=1

# Inform WEST where to find Python and our other scripts where to find WEST
export WEST_PYTHON=$(which python2.7)
if [[ -z "$WEST_ROOT" ]]; then
    echo "Must set environ variable WEST_ROOT"
    exit
fi

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    # The way we're calling this, it's $SLURM_SUBMIT_DIR, which is fine.
    export WEST_SIM_ROOT="$PWD"
fi

# Set simulation name
export SIM_NAME=$(basename $WEST_SIM_ROOT)
echo "simulation $SIM_NAME root is $WEST_SIM_ROOT"

# Export environment variables for the ZMQ work manager.
export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=300
export BASH=$SWROOT/bin/bash
export PERL=$SWROOT/usr/bin/perl
export ZSH=$SWROOT/bin/zsh
export IFCONFIG=$SWROOT/bin/ifconfig
export CUT=$SWROOT/usr/bin/cut
export TR=$SWROOT/usr/bin/tr
export LN=$SWROOT/bin/ln
export CP=$SWROOT/bin/cp
export RM=$SWROOT/bin/rm
export SED=$SWROOT/bin/sed
export CAT=$SWROOT/bin/cat
export HEAD=$SWROOT/bin/head
export TAR=$SWROOT/bin/tar
export AWK=$SWROOT/usr/bin/awk
export PASTE=$SWROOT/usr/bin/paste
export GREP=$SWROOT/bin/grep
export SORT=$SWROOT/usr/bin/sort
export UNIQ=$SWROOT/usr/bin/uniq
export HEAD=$SWROOT/usr/bin/head
export MKDIR=$SWROOT/bin/mkdir
export ECHO=$SWROOT/bin/echo
export DATE=$SWROOT/bin/date
export SANDER=$AMBERHOME/bin/sander
export PMEMD=$AMBERHOME/bin/pmemd.cuda
export CPPTRAJ=$AMBERHOME/bin/cpptraj
