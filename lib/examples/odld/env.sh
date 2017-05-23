if [[ -z "$WEST_ROOT" ]]; then
    echo "Must set environ variable WEST_ROOT"
    exit
fi

if [[ -z "$WEST_PYTHON" ]]; then
    export WEST_PYTHON=$(which python2.7)
fi
export WM_WORK_MANAGER=serial

