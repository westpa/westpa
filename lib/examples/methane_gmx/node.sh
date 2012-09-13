#!/bin/bash
unset MODULE_VERSION_STACK MODULE_VERSION MODULEPATH LOADEDMODULES MODULESHOME module
source /etc/profile.d/modules.sh

cd $WEMD_SIM_ROOT
source env.sh

set -x
$WEMD_ROOT/bin/w_run "$@"
