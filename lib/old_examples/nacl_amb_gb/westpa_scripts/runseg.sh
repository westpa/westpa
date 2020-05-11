#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

# Set up the run
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF
ln -sv $WEST_SIM_ROOT/amber_config/nacl.prm .

case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   parent segment
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/amber_config/md-continue.in > md.in
        ln -sv $WEST_PARENT_DATA_REF/seg.log ./parent.log
        ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory
        # $WEST_PARENT_DATA_REF contains the reference to the
        #   appropriate basis or initial state
        sed "s/RAND/$WEST_RAND16/g" \
          $WEST_SIM_ROOT/amber_config/md-genvel.in > md.in
        ln -sv $WEST_PARENT_DATA_REF ./parent.rst
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

# Propagate segment
$PMEMD -O -i md.in   -p nacl.prm  -c parent.rst \
          -r seg.rst -x seg.crd   -o seg.log    -inf seg.nfo

# Calculate progress coordinate
COMMAND="         trajin $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.crd\n"
COMMAND="$COMMAND distance na-cl :1@Cl- :2@Na+ out dist.dat\n"
COMMAND="$COMMAND go\n"
echo -e $COMMAND | $CPPTRAJ nacl.prm
cat dist.dat | tail -n +2 | awk '{print $2}' > $WEST_PCOORD_RETURN

# Output coordinates
if [ ${WEST_COORD_RETURN} ]; then
    COMMAND="         trajin  $WEST_CURRENT_SEG_DATA_REF/parent.rst\n"
    COMMAND="$COMMAND trajin  $WEST_CURRENT_SEG_DATA_REF/seg.crd\n"
    COMMAND="$COMMAND trajout $WEST_CURRENT_SEG_DATA_REF/seg.pdb\n"
    COMMAND="$COMMAND go\n"
    echo -e $COMMAND | $CPPTRAJ nacl.prm
    cat $WEST_CURRENT_SEG_DATA_REF/seg.pdb | grep 'ATOM' \
      | awk '{print $6, $7, $8}' > $WEST_COORD_RETURN
fi

# Output log
if [ ${WEST_LOG_RETURN} ]; then
    if [ $WEST_CURRENT_SEG_INITPOINT_TYPE == SEG_INITPOINT_CONTINUES ]; then
        cat $WEST_CURRENT_SEG_DATA_REF/parent.log \
          | awk '/RESULTS/ {p=1}; p; /A V E R A G E S/ {p=0}' | grep = \
          | tail -n 5 > $WEST_LOG_RETURN
    fi
    cat $WEST_CURRENT_SEG_DATA_REF/seg.log \
      | awk '/RESULTS/ {p=1}; p; /A V E R A G E S/ {p=0}' | grep = \
      >> $WEST_LOG_RETURN
fi

# Clean up
rm -f dist.dat md.in nacl.prm parent.log parent.rst seg.nfo seg.pdb
