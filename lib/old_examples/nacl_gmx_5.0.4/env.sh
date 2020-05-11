# This file defines where WEST and Gromacs can be found
# Modify to taste

# Set environment variables for Gromacs for convenience
export GROMPP=$(which grompp)
export MDRUN=$(which mdrun)
GMX_VERSION=$($MDRUN -version | grep 'Gromacs version' | awk '{print $4}' \
  | cut -c 1)
if [ $GMX_VERSION -eq 4 ]; then
    export G_DIST=$(which g_dist)
    export TRJCONV=$(which trjconv)
elif [ $GMX_VERSION -eq 5 ]; then
    export GMX=$(which gmx)
fi

# Inform WEST where to find Python and our other scripts where to find WEST
export WEST_PYTHON=$(which python2.7)
if [[ -z "$WEST_ROOT" ]]; then
    echo "Must set environ variable WEST_ROOT"
    exit
fi

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    export WEST_SIM_ROOT="$PWD"
fi

# Set simulation name
export SIM_NAME=$(basename $WEST_SIM_ROOT)
echo "simulation $SIM_NAME root is $WEST_SIM_ROOT"
