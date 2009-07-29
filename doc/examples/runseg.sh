function get_timing() {
    tail -50 $1 | perl -e \
'while($_ = <STDIN>) {
    if (m/Time:\s*(\S+)\s+(\S+)/) {
        print "$1";
        exit;
    } elsif (m/^real\s*(\d+)m([0-9\.]+)s/) {
        $secs = $1 * 60 + $2;
        print "$secs";
        exit;
    }
}'
}

source /opt/sci/gromacs/gromacs-4/bin/GMXRC.bash

if [ -z "$TMPDIR" ] ; then
    TMPDIR=/scratch/mzwier/$$
    mkdir -p $TMPDIR
fi

SRCDIR=/scratch/mzwier/we_test
WE_STATE=$SRCDIR/we.sqlite; export WE_STATE
PYTHONPATH="/home/mzwier/WE/workspace/we/src:$PYTHONPATH"; export PYTHONPATH
WE="python /home/mzwier/WE/workspace/we/src/we.py"; export WE

cd $TMPDIR
# What segment are we running?
eval `$WE nextseg`
WORKDIR=$TMPDIR/$WE_CURRENT_ITER/$WE_CURRENT_SEG_ID
mkdir -p $WORKDIR || exit 1
cd $WORKDIR

[ -z "$WE_CURRENT_ITER" ] && exit 1 

# Set up the run
cp $SIM_ROOT/{CH4.top,md.mdp,md.ndx} .
if [ -z "$WE_PARENT_SEG_DATA_REF" ] ; then
    # Bootstrap simulation
    ln -s $SIM_ROOT/eq2.gro ./startpos.gro
    ln -s $SIM_ROOT/eq2.trr ./parent.trr
else
    ln -s $WE_PARENT_SEG_DATA_REF/endpos.gro ./startpos.gro
    ln -s $WE_PARENT_SEG_DATA_REF/seg.trr ./parent.trr
fi

grompp -f md.mdp -c startpos.gro -t parent.trr -p CH4.top -o seg.tpr \
       &> setup.log || exit 1

(time mdrun -s seg.tpr -o seg.trr -x seg.xtc -c endpos.gro \
           -e seg.edr -g seg.log \
     ) &> mdrun.log || exit 1

# Get progress coordinate
echo "3 4" | g_dist -f seg.trr -s seg.tpr -n md.ndx &> g_dist.log
PCOORD=$(cat dist.xvg \
         | perl -e '@lines = <STDIN>; 
                    $lines[-1] =~ /\s*\S+\s+(\S+)/; 
                    print $1*(-10), "\n"')

# Save timing info
CPUTIME=$(get_timing mdrun.log)
$WE updateseg -c -t $CPUTIME -P $PCOORD || exit 1

mkdir -p $WE_CURRENT_SEG_DATA_REF
if [ -n "$WE_PARENT_SEG_DATA_REF" ] ; then
    ln -s $WE_PARENT_SEG_DATA_REF $WE_CURRENT_SEG_DATA_REF/parent
fi
cp seg.{tpr,trr,xtc,edr,log} endpos.gro *.log dist.xvg \
   $WE_CURRENT_SEG_DATA_REF

