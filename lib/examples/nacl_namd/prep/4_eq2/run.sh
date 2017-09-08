#!/bin/bash
cat ../../namd_config/toppar/toppar_water_ions.str \
  | grep -v -e "^set" \
  | grep -v -e "^if" \
  | grep -v -e "^WRNLEV"\
  | grep -v -e "^BOMLEV"\
  | grep -v -e "^SOD\s\s\s\sO"\
  > ../../namd_config/toppar/toppar_water_ions_for_namd.str
namd2 4_eq2.conf > 4_eq2.log

