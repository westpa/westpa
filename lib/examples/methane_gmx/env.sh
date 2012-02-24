# SIM_ROOT (conveniently?) defaults to the directory containing env.sh
source /etc/profile.d/modules.sh

if [[ -z "$WEMD_SIM_ROOT" ]]; then
    export WEMD_SIM_ROOT="$PWD"
fi
export SIM_NAME=$(basename $WEMD_SIM_ROOT)

echo "simulation $SIM_NAME root is $WEMD_SIM_ROOT"

case `hostname` in
    assign*)
        ulimit -c unlimited
        module load epd
        module unload gromacs
        module load gromacs/4.5
        #module load gromacs/4.5_git_static-gcc-4.5
        export WEMD_ROOT=$HOME/Eclipse/wemd_feat
        export WEMD_PYTHON=$(which python2.7)
    ;;

    ltc1*|compute8*|compute48*)
        module load epd
        module load gromacs/4.5
        export WEMD_ROOT=$HOME/wemd_feat
        export WEMD_PYTHON=$(which python2.7)
    ;;

    *)
        echo "unknown host in env.sh" > /dev/stderr
    ;;
esac
