#!/bin/bash

pmemd -O -i 5_run.in -o 5_run.out -c ../4_eq2/4_eq2.rst -p ../1_leap/nacl.parm7 -r 5_run.rst -ref ../4_eq2/4_eq2.rst

