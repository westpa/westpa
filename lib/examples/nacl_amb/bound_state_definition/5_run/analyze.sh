#!/bin/sh

COMMAND="         parm ../1_leap/nacl.parm7 \n"
COMMAND="$COMMAND trajin mdcrd \n"
COMMAND="$COMMAND distance na-cl :1@Na+ :2@Cl- out distance.dat\n"
COMMAND="$COMMAND go"

echo -e "$COMMAND" | cpptraj

python plot.py
