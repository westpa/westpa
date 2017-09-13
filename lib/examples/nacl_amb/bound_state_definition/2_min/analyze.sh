#!/bin/bash

CMD="     parm ../1_leap/nacl.parm7\n"
CMD="$CMD trajin 2_min.rst\n"
CMD="$CMD distance :1 :2 out dist.dat\n"
CMD="$CMD go \n"

echo -e "${CMD}" | cpptraj

cat dist.dat | tail -n 1 | awk '{print $2}'
