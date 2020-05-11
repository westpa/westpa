#!/bin/bash
pmemd -O -i 4_eq2.in -o 4_eq2.out -c ../3_eq1/3_eq1.rst -p ../1_leap/nacl.parm7 -r 4_eq2.rst -ref ../3_eq1/3_eq1.rst

