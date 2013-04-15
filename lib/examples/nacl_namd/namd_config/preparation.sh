#!/bin/bash

VMD=$(which vmd)
$VMD -dispdev text -e sys_build.tcl

rm -vr ../bstates
mkdir  ../bstates
cp     nacl.pdb ../bstates
outfile="../bstates/bstates.txt"
echo 0 1 nacl.pdb > $outfile
