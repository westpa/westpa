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
import numpy.ma as ma

from . import UncertMath
from . import BinCluster


def probAdjustEquil(binProb,rates,uncert,threshold=0.0,fullCalcClust=False,fullCalcBins=False):
    """This function adjusts bin pops in binProb using rates and uncert matrices
    fullCalcBins --> True for weighted avg, False for simple calc
    fullCalcClust --> True for weighted avg, False for simple calc
    threshold --> minimum weight (relative to max) for another value to be averaged
            only matters if fullCalcBins == True (or later perhaps if fullCalcClust == True)
    """
    
    
    # Check that rate matrix is square
    Ni,Nj = rates.shape
    if Ni != Nj:
        print('\nWARNING: Not a square matrix!\n')

    zi = np.where(binProb == 0.0)[0]  # indices of bins with zero probability
    
    rates_uncert = UncertMath.UncertContainer(rates,rates - uncert,rates + uncert)
    
    # STEP 1a: Create matrix of ratios of probabilities based on DIRECT estimates
    # that is, ij element is p_i / p_j = k_ji / k_ij
    
    ratios_direct = rates_uncert.transpose() / rates_uncert  

    # STEP 1b: Create averaged matrix of ratios of probabilities based on both direct and indirect estimates
    # Indirect means '3rd bin' estimates: p_i / p_j = ( k_ki / k_ik ) ( k_jk / k_kj )
    # Turns out this is not helpful, so generally set fullCalcBins = 0 
    if fullCalcBins:
        # Calculate indirect ratios using Einstein Summation convention where
        # ratios_indirect_kij  = ( k_ki / k_ik ) ( k_jk / k_kj ) = ratios_direct_ik * ratios_direct_kj
        ri_vals = np.einsum('ik,kj->kij',ratios_direct.vals,ratios_direct.vals)
        ri_min = np.einsum('ik,kj->kij',ratios_direct.dmin,ratios_direct.dmin)
        ri_max = np.einsum('ik,kj->kij',ratios_direct.dmax,ratios_direct.dmax)
        ratios_indirect = UncertMath.UncertContainer(ri_vals,ri_min,ri_max,mask=ratios_direct.vals.mask)

        # Threshold indirect ratios 
        ti = ratios_indirect.wt < ratios_direct * threshold
        ratios_indirect.mask = ti
        ratios_indirect.update_mask()

        ratios_indirect.concatenate(ratios_direct,axis=0) 
        ratios_average = ratios_indirect.weighted_average(axis=0)
 
    else:
        ratios_average = ratios_direct.weighted_average(axis=0,expaxis=0)
    

    # STEP 2: Form clusters

    # STEP 2a: Sort probability ratios based on uncertainty
    # Sort uncertainties of ratios_average subject to the convention that p_i < p_j
    
    i,j = np.triu_indices(Ni,1) # indices of ij pairs where i != j

    # Remove pairs that include a bin that has zero probability
    nzi = (binProb[i] != 0.0) & (binProb[j] != 0.0)
    i = i[nzi]
    j = j[nzi]

    vals = ma.vstack((ratios_average.vals[i,j],ratios_average.vals[j,i]))
    ias = ma.argsort(vals,axis=0,fill_value=np.inf)
    
    ordered_ind = np.vstack((i,j))
    flip_ind = np.nonzero(ias[0,:]) # Find pairs in which to select ji rather than ij
    ordered_ind[:,flip_ind[0]] = ordered_ind[:,flip_ind[0]][::-1]
    
    iind = ordered_ind[0,:]
    jind = ordered_ind[1,:]
    uncertij = ratios_average.uncert[iind,jind] # Get the uncert for ij pairs

    count = uncertij.count() # Count of the unmasked uncertainties
    ias = ma.argsort(uncertij,fill_value=np.inf) # Get the indices that would sort uncertij
    iind = iind[ias[:count]] # Sort the indices excluding masked/undefined values
    jind = jind[ias[:count]]


    # STEP 2b: Create ClusterList object and cluster bins
    clusters = BinCluster.ClusterList(ratios_average,Ni)

    if fullCalcClust:
        clusters.join((iind,jind))
    else:
        clusters.join_simple((iind,jind))

    total_prob = 0.0  # total probability in all clusters
    for cid in clusters.cluster_contents:
        binlist = list(clusters.cluster_contents[cid])
        if len(binlist):
            prob_cluster = binProb[binlist].sum()
            total_prob += prob_cluster

            binProb[binlist] = prob_cluster * clusters.bin_data[binlist].vals

    binProb[zi] = 0.0 # re-zero bins that previously had zero prob
    #for bi,p in enumerate(binProb):
    #    print('bin: {} -- {}'.format(bi,p))
    print('.........Total Probability: {}'.format(binProb.sum()))
            
    
