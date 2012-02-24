#!/bin/bash
unset MODULE_VERSION_STACK MODULE_VERSION MODULEPATH LOADEDMODULES MODULESHOME module
source /etc/profile.d/modules.sh

cd $SIM_ROOT
source env.sh

$WEMD_ROOT/bin/w_run "$@"
