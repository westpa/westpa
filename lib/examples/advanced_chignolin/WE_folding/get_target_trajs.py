import os
import sys
import re
import subprocess
import h5py
import numpy as np


MinIter = int(sys.argv[1])
MaxIter = int(sys.argv[2])
TargetBinVal = 0.5  ## This is the folding target state RMSD


DataIn = h5py.File("./west.h5", 'r')

FileOut = open("target_trajs.dat", 'w')
print>>FileOut, "iteration ID\t segment ID\n"


for i in range(MinIter,MaxIter+1) :
   if i%100 == 0 : print "Iter ", i
   ITER = "{0:08}".format(i)
   SegNum = int(DataIn["/iterations/iter_"+ITER+"/pcoord"].shape[0])

   for j in range(SegNum) :
      if float(DataIn["/iterations/iter_"+ITER+"/pcoord"].value[j][1]) < TargetBinVal :
         print>>FileOut, i, "\t", j 

FileOut.close()
