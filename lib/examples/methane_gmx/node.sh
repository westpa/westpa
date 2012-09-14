#!/bin/bash
unset MODULE_VERSION_STACK MODULE_VERSION MODULEPATH LOADEDMODULES MODULESHOME module
source /etc/profile.d/modules.sh

cd $WEST_SIM_ROOT
source env.sh

set -x
$WEST_ROOT/bin/w_run "$@"
