# This file defines where WEST and Amber can be found.
# Modify to taste

# Inform WEST where to find Amber.
# If AMBERHOME isn't set, it's unlikely that the libraries can be found, regardless of whether cpptraj/pmemd exist in the path.
if [[ -z "$AMBERHOME" ]]; then
    echo "AMBERHOME environ variable not set.  Source amber.sh, or modify env.sh such that WESTPA knows where AMBER is."
    exit
fi

export PATH=$AMBERHOME/bin:$PATH

# Set environment variables for Amber for convenience
export PMEMD=$(which pmemd)
export CPPTRAJ=$(which cpptraj)

# Inform WEST where to find Python and our other scripts where to find WEST
export WEST_PYTHON=$(which python3)

if [[ -z "$WEST_ROOT" ]]; then
    echo "Please set the environment variable WEST_ROOT"
    exit
fi

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    export WEST_SIM_ROOT="$PWD"
fi

# Set simulation name
export SIM_NAME=$(basename $WEST_SIM_ROOT)
echo "simulation $SIM_NAME root is $WEST_SIM_ROOT"
