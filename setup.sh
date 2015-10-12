#!/bin/bash

function checkout_remote() {
    directory="$1"; shift
    repo="$1"; shift
    revision="$1"; shift

    pushd lib    
    if ! [[ -d $directory ]] ; then

        git clone "$repo" || exit 1
    fi
    pushd $directory
      git pull origin master || exit 1
      git checkout "$revision" || exit 1
    popd
    popd
}

if [[ -z "$WEST_PYTHON" ]] ;  then
    WEST_PYTHON=python2.7
fi

find . -name \*.so -print0 | xargs -0 rm &> /dev/null

#checkout_remote blessings  git://github.com/erikrose/blessings.git d3ba51c5870d599b40b387ac6703805c3e23d292 || exit 1
#checkout_remote west_tools git://github.com/westpa/west_tools.git  15ac2534bb94c7386938192bd6d70dbf8401bf24 || exit 1
#checkout_remote wwmgr      git://github.com/westpa/wwmgr.git       06199824dbadd5e6ea184ac1f849339dff05e652 || exit 1

if [[ -d lib/h5py ]] ; then
    echo "using custom h5py located in $PWD/lib/h5py"

    if [[ -z "$HDF5" ]] ; then
        export HDF5=$(ldd `which h5ls` | fgrep libhdf5 | awk '{print $3;}' | sed 's~/lib.*~~')
    fi
    echo "using HDF5 from $HDF5"
    echo "(if this is not what you expect, set HDF5=... prior to running this script)"

    pushd lib/h5py
      cd h5py
        $WEST_PYTHON api_gen.py
      cd ..
      $WEST_PYTHON setup.py build_ext --inplace --hdf5=$HDF5
    popd
fi

for d in lib/west_tools; do
    cd $d
    $WEST_PYTHON setup.py build_ext --inplace
    cd -
done

