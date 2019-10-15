#!/bin/bash
source env.sh

SFX=.d$$
mv traj_segs{,$SFX}
mv seg_logs{,$SFX}
mv istates{,$SFX}
rm -Rf traj_segs$SFX seg_logs$SFX istates$SFX & disown %1
rm -f system.h5 west.h5 seg_logs.tar
mkdir seg_logs traj_segs istates

# Create our reference state.
echo -e "4 \n" | $GMX trjconv -f bstates/0/seg.gro -s bstates/0/seg.tpr -n 18-crown-6-K+.ndx -o initial_bin_struct.pdb || exit 1
echo "2 9" | $GMX trjconv -fit rot+trans -s bound_state.tpr -f initial_bin_struct.pdb -o initial_bin_struct_align.pdb || exit 1

cat initial_bin_struct_align.pdb | grep '^ATOM' | grep K\+ > 18-crown-6-K+.pdb
rm initial_bin_struct.pdb initial_bin_struct_align.pdb

BSTATE_ARGS="--bstates-from BASIS_STATES.single"

$WEST_ROOT/bin/w_init $BSTATE_ARGS $TSTATE_ARGS --segs-per-state 1 --serial "$@"
