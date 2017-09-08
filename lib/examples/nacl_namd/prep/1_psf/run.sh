#!/bin/bash
cat ../../namd_config/toppar/toppar_water_ions.str \
  | grep -v -e "^set" \
  | grep -v -e "^if" \
  | grep -v -e "^WRNLEV"\
  | grep -v -e "^BOMLEV"\
  | grep -v -e "^SOD\s\s\s\sO"\
  > ../../namd_config/toppar/toppar_water_ions_for_namd.str

vmd -dispdev text -e prep.pgn

