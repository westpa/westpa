#!/bin/bash
pmemd -O -i 3_eq1.in -o 3_eq1.out -c ../2_min/2_min.rst -p ../1_leap/nacl.parm7 -r 3_eq1.rst -ref ../2_min/2_min.rst

