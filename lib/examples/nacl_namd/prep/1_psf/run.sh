#!/bin/bash

# Make sure the CHARMM force field is in the right spot:
if [ ! -d ../../namd_config/toppar ]; then
  echo "CHARMM force field not found! Download the CHARMM36 force field, untar "
  echo "it, and place the 'toppar' directory in ../../namd_config."
  exit 1
fi

cat ../../namd_config/toppar/toppar_water_ions.str \
  | grep -v -e "^set" \
  | grep -v -e "^if" \
  | grep -v -e "^WRNLEV"\
  | grep -v -e "^BOMLEV"\
  | grep -v -e "^SOD\s\s\s\sO"\
  > ../../namd_config/toppar/toppar_water_ions_for_namd.str

vmd -dispdev text -e prep.pgn

