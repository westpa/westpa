#! /bin/bash -l

if [ -z "$WEST_ROOT" ]; then
    echo "Please define a $WEST_ROOT"
    read WEST_ROOT
fi

chmod +x w_ipython
cp w_ipython $WEST_ROOT/bin
cp w_ipython.py $WEST_ROOT/lib/west_tools

echo "Done!"