module load epd &> /dev/null

export WEMD_ROOT=$(readlink -f $PWD/../../..)
export WEMD_PYTHONPATH=$HOME/Eclipse/wwmgr
export WEMD_PYTHON=$EPD_ROOT/bin/python2.7
export WWMGR_WORK_MANAGER=serial


