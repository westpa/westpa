#!/bin/bash
set -x
cd $WEST_SIM_ROOT/bstates || exit 1
cat pcoord.init > $WEST_PCOORD_RETURN
