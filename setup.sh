#! /bin/bash -l

if [ -z "$WEST_ROOT" ]; then
    echo "Please define a $WEST_ROOT"
    read WEST_ROOT
fi

chmod +x w_ipython
cp w_ipython.py $WEST_ROOT/lib/west_tools
ln -sv $WEST_ROOT/bin/west $WEST_ROOT/bin/w_ipython

echo "Done!"
