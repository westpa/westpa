module load epd &> /dev/null

export WEST_ROOT=$(readlink -f $PWD/../../..)
export WEST_PYTHONPATH=$HOME/Eclipse/wwmgr
export WEST_PYTHON=$EPD_ROOT/bin/python2.7
export WWMGR_WORK_MANAGER=serial


