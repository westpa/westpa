#!/bin/bash

$WEST_ROOT/bin/w_assign --states-from-function system.gen_state_labels
$WEST_ROOT/bin/w_pdist -o pdist.h5
$WEST_ROOT/bin/plothist average --first-iter 1000 --log10 pdist.h5

$WEST_ROOT/bin/w_kinetics trace
$WEST_ROOT/bin/w_kinavg trace -e cumulative --window-frac .5 --step-iter 10
$WEST_ROOT/bin/ploterr kinavg

