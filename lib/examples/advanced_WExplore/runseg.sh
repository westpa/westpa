#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT
source env.sh

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd $WEST_CURRENT_SEG_DATA_REF || exit 1

if [[ "$USE_LOCAL_SCRATCH" == "1" ]] ; then
    # make scratch directory
    WORKDIR=/$SCRATCH/$USER/$WEST_JOBID/$WEST_CURRENT_SEG_DATA_REF
    $SWROOT/bin/mkdir -p $WORKDIR || exit 1
    cd $WORKDIR || exit 1
    STAGEIN="$SWROOT/bin/cp -avL"
else
    STAGEIN="$SWROOT/bin/ln -sv"
fi

function cleanup() {
    # Clean up
    if [[ "$USE_LOCAL_SCRATCH" == "1" ]] ; then
        $SWROOT/bin/cp *.{trr,edr,gro,pdb} $WEST_CURRENT_SEG_DATA_REF || exit 1
        cd $WEST_CURRENT_SEG_DATA_REF
        #$SWROOT/bin/rm -Rf $WORKDIR
    else
        $SWROOT/bin/rm -f *.itp *.mdp *.ndx *.top state.cpt
        $SWROOT/bin/rm -f none.xtc whole.xtc
    fi
}

trap cleanup EXIT

# Set up the run
$STAGEIN $WEST_SIM_ROOT/{$TOP,$NDX,*.mdp,*.itp,$REF} .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        $STAGEIN $WEST_PARENT_DATA_REF/seg.gro ./parent.gro
        $STAGEIN $WEST_PARENT_DATA_REF/seg.trr ./parent.trr
        $STAGEIN $WEST_PARENT_DATA_REF/seg.edr ./parent.edr
        $STAGEIN $WEST_PARENT_DATA_REF/imaged_ref.gro ./parent_imaged.gro
        $GMX grompp -f md.mdp -c parent.gro -t parent.trr -p $TOP -o seg.tpr -n $NDX -maxwarn 2 -e parent.edr || exit 1
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
        $STAGEIN $WEST_PARENT_DATA_REF.gro ./initial.gro
        $STAGEIN $WEST_PARENT_DATA_REF.trr ./initial.trr
        $STAGEIN $WEST_PARENT_DATA_REF.edr ./initial.edr
        $STAGEIN $WEST_PARENT_DATA_REF.gro ./parent_imaged.gro
        $GMX grompp -f md.mdp -c initial.gro -p $TOP -o seg.tpr -n $NDX -maxwarn 2 -t initial.trr -e initial.edr|| exit 1
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac
    
# Propagate segment
export OMP_NUM_THREADS=1
$GMX mdrun -s seg.tpr -o seg.trr -x seg.xtc -c seg.gro \
      -e seg.edr -g seg.log -nt 1 \
      || exit 1

# Image trajectory (to correct for jumps out of the box).
echo -e "0 \n" | $GMX trjconv    -f seg.trr  -s parent_imaged.gro  -n $NDX -o none.xtc  -pbc none || exit 1
echo -e "0 \n" | $GMX trjconv    -f none.xtc  -s parent_imaged.gro  -n $NDX -o whole.xtc -pbc whole || exit 1
echo -e "0 \n" | $GMX trjconv    -f whole.xtc  -s parent_imaged.gro  -n $NDX -o nojump.xtc -pbc nojump || exit 1
# This next command needs the exact timepoint for -b, in ps.  As this uses a 2 ps tau, to get the final frame, we specify 2.
echo -e "0 \n" | $GMX trjconv    -f nojump.xtc  -s seg.tpr -n $NDX -o imaged_ref.gro -b 2 || exit 1
echo -e "4 \n" | $GMX trjconv    -f nojump.xtc  -s seg.tpr -n $NDX -o pcoord.pdb || exit 1

# Copy in the imaged trajectory as the progress coordinate.  We'll use a python
# pcoord loader to analyze the RMSD and go from there.
echo "2 9" | $GMX trjconv -fit rot+trans -s bound_state.tpr -f pcoord.pdb -o pcoord_align.pdb
cat pcoord_align.pdb | grep '^ATOM' | grep K\+ > $WEST_PCOORD_RETURN || exit 1
