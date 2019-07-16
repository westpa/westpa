# Copyright (C) 2013 Joshua L. Adelman, Carsen A. Stringer and Daniel M. Zuckerman
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np


from . import UncertMath

class ClusterList(object):
    def __init__(self,ratios,nbins):
        super(ClusterList, self).__init__()
        self.nbins = nbins
        self.ratios = ratios

        # Create an array to hold bin assignments and initial set all to -1 == not clustered
        self.bin_assign = np.empty((nbins,),dtype=np.int) 
        self.bin_assign.fill(-1)                         

        # Initialize an ucert container to hold per bin information; initially mask all elements
        # note: masking is implicit since rates set to 0
        dum_data = np.zeros((nbins,))
        self.bin_data = UncertMath.UncertContainer(dum_data.copy(),dum_data.copy(),dum_data.copy())

        self.cluster_id = 0                     # ID of newest cluster
        self.cluster_contents = {}              # Dictionary containing sets of bin ids with keys = cluster ids
    
    def join(self,pairs):
        """ Join clusters given a tuple (i,j) of bin pairs
        """

        for i,j in zip(*pairs):
            # Both bins not joined
            if self.bin_assign[i] == -1 and self.bin_assign[j] == -1:
                # Create new cluster
                self.bin_assign[i] = self.cluster_id
                self.bin_assign[j] = self.cluster_id
                
                self.cluster_contents[self.cluster_id] = {i,j}
                self.cluster_id += 1

                rij = self.ratios[i,j]
                denom = rij + 1.0
                self.bin_data[i] = rij/denom          # relative probability for bin i
                self.bin_data[j] = denom.recip()      # relative probability for bin j

            # Only one bin previously assigned to a cluster
            elif  self.bin_assign[i] == -1 or self.bin_assign[j] == -1:
                if self.bin_assign[i] == -1:
                    idum,jdum = i,j
                else:
                    idum,jdum = j,i
                
                jclust = self.bin_assign[jdum]
                jclust_mid = np.where(self.bin_assign == jclust)[0]    # index of bins in jclust

                rik = self.ratios[idum,jclust_mid]
                pk = self.bin_data[jclust_mid]
                piTmp = rik * pk     # estimate for p_idum / P_cluster based on 'path' through bin k
                                     # Note that here P_cluster is value before addition of bin idum
                piTmp_avg = piTmp.weighted_average(axis=0)

                # now, compute relative prob of each bin in *new* cluster (including bin idum)   
                denom = piTmp_avg + 1.0
                self.bin_data[idum] = piTmp_avg / denom
                
                # Update bins already in cluster
                self.bin_data[jclust_mid] = self.bin_data[jclust_mid] / denom

                # Move bin idum into cluster jclust
                self.bin_assign[idum] = jclust
                self.cluster_contents[jclust].update({idum})

            # Both bins previously assigned to different cluster; Join clusters
            elif not self.bin_assign[i] == self.bin_assign[j]:
                iclust = self.bin_assign[i]
                jclust = self.bin_assign[j]

                iclust_mid = np.where(self.bin_assign == iclust)[0]  # indx of bins in cluster i
                jclust_mid = np.where(self.bin_assign == jclust)[0]  # indx of bins in cluster j
                
                niclust = iclust_mid.size
                njclust = jclust_mid.size
                
                dum_data = np.zeros((niclust*njclust,))
                ij_cluster_ratio = UncertMath.UncertContainer(dum_data.copy(),dum_data.copy(),dum_data.copy())
                 
                for count,im in enumerate(iclust_mid): 
                    rij = self.ratios[im,jclust_mid]
                    pi = self.bin_data[im]
                    pj = self.bin_data[jclust_mid]

                    ij_cluster_ratio[count*njclust:(count+1)*njclust] = rij * pj / pi

                ij_cluster_ratio = ij_cluster_ratio.weighted_average(axis=0)

                idenom = ij_cluster_ratio.recip() + 1.0
                jdenom = ij_cluster_ratio + 1.0
                
                self.bin_data[iclust_mid] = self.bin_data[iclust_mid] / idenom
                self.bin_data[jclust_mid] = self.bin_data[jclust_mid] / jdenom
                
                # Join all bins in cluster j into cluster iclust
                self.bin_assign[jclust_mid] = iclust

                # Move contents of jclust into iclust
                self.cluster_contents[iclust].update(self.cluster_contents[jclust])
                
                # Clear contents of jclust
                self.cluster_contents[jclust].clear()

                if len(self.cluster_contents[iclust]) == self.nbins:
                    break
            
    def join_simple(self,pairs):
        """ Join clusters using direct ratios given a tuple (i,j) of bin pairs
        """

        for i,j in zip(*pairs):
            # Both bins not joined
            if self.bin_assign[i] == -1 and self.bin_assign[j] == -1:
                # Create new cluster
                self.bin_assign[i] = self.cluster_id
                self.bin_assign[j] = self.cluster_id
                
                self.cluster_contents[self.cluster_id] = {i,j}
                self.cluster_id += 1

                rij = self.ratios[i,j]
                denom = rij + 1.0
                self.bin_data[i] = rij/denom          # relative probability for bin i
                self.bin_data[j] = denom.recip()      # relative probability for bin j

            # Only one bin previously assigned to a cluster
            elif  self.bin_assign[i] == -1 or self.bin_assign[j] == -1:
                if self.bin_assign[i] == -1:
                    idum,jdum = i,j
                else:
                    idum,jdum = j,i
                
                jclust = self.bin_assign[jdum]
                rik = self.ratios[idum,jdum]
                pk = self.bin_data[jdum]
                piTmp = rik * pk     # estimate for p_idum / P_cluster based on 'path' through bin k
                                     # Note that here P_cluster is value before addition of bin idum
                # now, compute relative prob of each bin in *new* cluster (including bin idum)   
                denom = piTmp + 1.0
                self.bin_data[idum] = piTmp / denom
                
                # Update bins already in cluster
                jclust_mid = np.where(self.bin_assign == jclust)    # index of bins in jclust
                self.bin_data[jclust_mid] = self.bin_data[jclust_mid] / denom

                # Move bin idum into cluster jclust
                self.bin_assign[idum] = jclust
                self.cluster_contents[jclust].update({idum})

            # Both bins previously assigned to different cluster; Join clusters
            elif not self.bin_assign[i] == self.bin_assign[j]:
                iclust = self.bin_assign[i]
                jclust = self.bin_assign[j]
                rij = self.ratios[i,j]
                pi = self.bin_data[i]
                pj = self.bin_data[j]
                ij_cluster_ratio = rij * pj / pi
                idenom = ij_cluster_ratio.recip() + 1.0
                jdenom = ij_cluster_ratio + 1.0

                iclust_mid = np.where(self.bin_assign == iclust)
                self.bin_data[iclust_mid] = self.bin_data[iclust_mid] / idenom

                jclust_mid = np.where(self.bin_assign == jclust)
                self.bin_data[jclust_mid] = self.bin_data[jclust_mid] / jdenom
                
                # Join all bins in cluster j into cluster iclust
                self.bin_assign[jclust_mid] = iclust

                # Move contents of jclust into iclust
                self.cluster_contents[iclust].update(self.cluster_contents[jclust])
                
                # Clear contents of jclust
                self.cluster_contents[jclust].clear()

                if len(self.cluster_contents[iclust]) == self.nbins:
                    break
            
        
