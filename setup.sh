#!/bin/bash

for d in src lib/west_tools; do
    cd $d
    find . -name \*.so -print0 | xargs -0 rm
    python setup.py build_ext --inplace
    cd -
done

