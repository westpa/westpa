# This file defines where WEST and NAMD can be found
# Modify to taste

# Inform WEST where to find Python and our other scripts where to find WEST
export WEST_PYTHON=$(which python2.7)
export WEST_ROOT=$(readlink -f $PWD/../../..)

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    export WEST_SIM_ROOT="$PWD"
fi
# Set environment variable for NAMD for convenience
export NAMD=$(which namd2)

# Set simulation name
export SIM_NAME=$(basename $WEST_SIM_ROOT)
echo "simulation $SIM_NAME root is $WEST_SIM_ROOT"
