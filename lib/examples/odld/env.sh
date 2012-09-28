source /etc/profile.d/modules.sh
module load epd &> /dev/null

export WEST_ROOT=$(readlink -f $PWD/../../..)
export WEST_PYTHON=$EPD_ROOT/bin/python2.7
export WM_WORK_MANAGER=serial

