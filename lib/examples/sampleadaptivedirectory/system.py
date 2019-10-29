from __future__ import print_function, division
import numpy
from west.propagators import WESTPropagator
from west.systems import WESTSystem
from westpa.binning import RectilinearBinMapper
from westpa.binning import FuncBinMapper
from westpa.binning import RecursiveBinMapper

import logging
log = logging.getLogger('westpa.rc')

PI = numpy.pi
from numpy import *
pcoord_dtype = numpy.float32
#THESE ARE THE THINGS YOU SHOULD CHANGE
bintargetcount=4 #number of walkers per bin
lengthpcoord=101 #the length of your progress coordinate 
numberofdim=2  # number of dimensions
binsperdim=8   # excludes the minimum and maximum case. You will have binsperbin^numberofdim+numberofdim*2 bins total
PCA=True       # choose to do principal component analysis
maxcap=[inf,inf]	#these are in the order of the dimensions left is first dimension and right is second dimension
mincap=[-inf,-inf]
targetstate=[None, None]    #enter boundaries for target state or None if there is no target state in that dimension
targetstatedirection=-1  #if your target state is meant to be greater that the starting pcoor use 1 or else use -1
activetarget=0		#if no target state make this zero
def function_map(coords, mask, output):
    	varcoords=copy(coords)
	originalcoords=copy(coords)
    	if PCA and len(output)>1:
    	    colavg=mean(coords, axis=0)
    	    for i in range(len(coords)):
    	    	for j in range(len(coords[i])):
    	    		varcoords[i][j]=coords[i][j]-colavg[j]
    	    covcoords=cov(transpose(varcoords))
    	    eigval, eigvec=linalg.eigh(covcoords)
	    eigvec=eigvec[:,argmax(absolute(eigvec),axis=1)]
	    for i in range(len(eigvec)):
		if eigvec[i,i]<0:
			eigvec[:,i]=-1*eigvec[:,i]
    	    for i in range(numberofdim):
    	    	for j in range(len(output)):
    	    		coords[j][i]=dot(varcoords[j],eigvec[:,i])
    	maxlist=[]
    	minlist=[]
    	for n in range(numberofdim):
    	    maxlist.append(amax(coords[:,n]))
    	    minlist.append(amin(coords[:,n]))
    	for i in range(len(output)):
    	    holder=2*numberofdim
    	    for n in range(numberofdim):
		    if (activetarget==1) and targetstate[n] is not None:
		    	if (originalcoords[i,n]*targetstatedirection) >= (targetstate[n]*targetstatedirection):
			    	holder=binsperdim**numberofdim+numberofdim*2
		    if (holder==binsperdim**numberofdim+numberofdim*2):
			    n=numberofdim
    	            elif coords[i,n]>=maxlist[n] or originalcoords[i,n]>=maxcap[n]:
    	                    holder= 2*n
    	                    n=numberofdim
    	            elif coords[i,n]<=minlist[n] or originalcoords[i,n]<=mincap[n]:
    	                    holder =2*n+1
    	                    n=numberofdim
    	    if holder==2*numberofdim:
    	    	for j in range(numberofdim):
    	    		holder = holder + (digitize(coords[i][j],linspace(minlist[j],maxlist[j],binsperdim+1))-1)*(binsperdim)**j
    	    output[i]=holder
    	return output

class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim = numberofdim
        self.pcoord_len = 101
	self.pcoord_dtype = numpy.float32 
        self.bin_mapper = FuncBinMapper(function_map, binsperdim**numberofdim+numberofdim*2+activetarget) #Changed binsperbin to binsperdim 
        self.bin_target_counts = numpy.empty((self.bin_mapper.nbins,), numpy.int_)
        self.bin_target_counts[...] = bintargetcount
