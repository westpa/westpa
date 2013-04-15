# sys_build.tcl
#   Generates psf file for NaCl
########################################### MODULES, SETTINGS, AND DEFAULTS ############################################
package     require     psfgen
################################################## GENERATE TOPOLOGY ###################################################
resetpsf
topology    toppar_water_ions.str
segment     1           {pdb nacl.pdb}
coordpdb    nacl.pdb    1
writepsf    nacl.psf
exit
