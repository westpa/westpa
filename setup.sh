#!/bin/bash


if [[ -d lib/h5py ]] ; then
    echo "using custom h5py located in $PWD/lib/h5py"

    if [[ -z "$HDF5" ]] ; then
        export HDF5=$(ldd `which h5ls` | fgrep libhdf5 | awk '{print $3;}' | sed 's~/lib.*~~')
    fi
    echo "using HDF5 from $HDF5"
    echo "(if this is not what you expect, set HDF5=... prior to running this script)"

    pushd lib/h5py
      cd h5py
        python api_gen.py
      cd ..
      python setup.py build_ext --inplace --hdf5=$HDF5
    popd
fi

for d in src lib/west_tools; do
    cd $d
    find . -name \*.so -print0 | xargs -0 rm
    python setup.py build_ext --inplace
    cd -
done

