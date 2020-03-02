#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

function cleanup() {
    rm -f pcoord.$$.pdb
    rm -f pcoord_align.$$.pdb
}

trap cleanup EXIT

# Get progress coordinate

echo $NDX
echo -e "4 \n" | $GMX trjconv -f $WEST_STRUCT_DATA_REF.gro  -s $WEST_STRUCT_DATA_REF.tpr -n $NDX -o pcoord.$$.pdb || exit 1
echo "2 9" | $GMX trjconv -fit rot+trans -s bound_state.tpr -f pcoord.$$.pdb -o pcoord_align.$$.pdb
cat pcoord_align.$$.pdb | grep '^ATOM' | grep K\+ > $WEST_PCOORD_RETURN

if [ -n "$SEG_DEBUG" ] ; then
    head -v $WEST_PCOORD_RETURN
fi


