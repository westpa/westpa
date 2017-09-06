#!/bin/sh

# Check for gmx
GMX=$(which gmx)
if [ ! -n "${GMX}" ]; then
  echo "Gromacs not found! Make sure your Gromacs installation is in your path."
  exit 1
fi
GROMACS_BIN=$(dirname ${GMX})
GROMACS_DIR=${GROMACS_BIN%/bin}
echo "Using force field files from Gromacs installation at $GROMACS_DIR"

#aminoacids.rtp
cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/aminoacids.rtp \
  | sed -n '/\[\sbondedtypes\s\]/,/^$/p'\
  > tip3p_ionsjc2008.ff/aminoacids.rtp

cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/aminoacids.rtp \
  | sed -n '/\[\sHOH\s\]/,/^$/p'\
  >> tip3p_ionsjc2008.ff/aminoacids.rtp

cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/aminoacids.rtp \
  | sed -n '/\[\sCL\s\]/,/^$/p'\
  >> tip3p_ionsjc2008.ff/aminoacids.rtp

cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/aminoacids.rtp \
  | sed -n '/\[\sNA\s\]/,/^$/p'\
  >> tip3p_ionsjc2008.ff/aminoacids.rtp

#atomtypes.atp
cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/atomtypes.atp \
  | grep -e "HW\s\|OW\s\|Cl\s\|Na\s" \
  > tip3p_ionsjc2008.ff/atomtypes.atp

#ffbonded.itp
echo "[ bondtypes ]" > tip3p_ionsjc2008.ff/ffbonded.itp
echo "; i    j  func       b0          kb" >> tip3p_ionsjc2008.ff/ffbonded.itp
cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/ffbonded.itp \
  | grep -e '\s\+OW\sHW\|\s\+HW\s\HW' \
  >> tip3p_ionsjc2008.ff/ffbonded.itp

echo "" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo "[ constrainttypes ]" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo "" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo "[ angletypes ]" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo ";  i    j    k  func       th0       cth" >> tip3p_ionsjc2008.ff/ffbonded.itp
cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/ffbonded.itp \
  | grep -e 'HW\s\sOW\s\sHW' \
  >> tip3p_ionsjc2008.ff/ffbonded.itp
echo "" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo "[ dihedraltypes ]" >> tip3p_ionsjc2008.ff/ffbonded.itp
echo ";i  j   k  l	 func      phase      kd      pn" >> tip3p_ionsjc2008.ff/ffbonded.itp

#ffnonbonded.itp
echo "[ atomtypes ]" > tip3p_ionsjc2008.ff/ffnonbonded.itp
echo "; name      at.num  mass     charge ptype  sigma      epsilon" >> tip3p_ionsjc2008.ff/ffnonbonded.itp
cat ${GROMACS_DIR}/share/gromacs/top/amber03.ff/ffnonbonded.itp\
  | grep -e "HW\s\|OW\s" \
  >> tip3p_ionsjc2008.ff/ffnonbonded.itp

echo "Cl          17      35.45    0.0000  A   4.47766e-01  1.48913e-01" >> tip3p_ionsjc2008.ff/ffnonbonded.itp
echo "Na          11      22.99    0.0000  A   2.43928e-01  3.65846e-01" >> tip3p_ionsjc2008.ff/ffnonbonded.itp 

#forcefield.itp
# already in directory

#ions.itp
# already in directory

#tip3p.itp
cp ${GROMACS_DIR}/share/gromacs/top/amber03.ff/tip3p.itp tip3p_ionsjc2008.ff/

#watermodels.dat
# already in directory
