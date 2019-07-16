# Copyright (C) 2013 Joshua L. Adelman
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
import itertools
import scipy.optimize


def solve_steady_state(T, U, target_bins_index):
        ntarget = len(target_bins_index)          # Number of target states
        nstates = T.shape[0]                    # Number of total states

        # Number of active states
        nsolve = nstates - ntarget

        # list of active states
        sactive = sorted(list(set(range(nstates)) - set(target_bins_index)))

        W = np.zeros((nsolve,nsolve))
        W_unc = np.zeros((nsolve,nsolve))
        S = np.zeros((nsolve,))

        for iiw,iit in zip(range(nsolve),sactive):
            for jjw,jjt in zip(range(nsolve),sactive):
                W[iiw,jjw] = T[jjt,iit]
                W_unc[iiw,jjw] = U[jjt,iit]

        for ii in range(nsolve):
            W[ii,ii] = 0.0
            for jj in range(nstates): #nstates
                if jj != ii:
                    W[ii,ii] -= T[ii,jj]

        #we have n equations, but only n-1 linearly independent equations
        #set one equation to constraint sum(Pi)=1

        # Pick row that contains the largest uncertainty to replace
        ii = np.unravel_index(np.argmax(W_unc),W_unc.shape)[1]

        S[ii] = 1.0
        W[ii,:] = 1.0

        #P = np.linalg.solve(W,S)
        try:
            P,err = scipy.optimize.nnls(W,S)
        except RuntimeError:
            print('Solve did not converge')
            return None

        # There are some instances where a single bin has its prob set 
        # to 1.0; In this case do not reweight
        if np.allclose(np.max(P), 1.0):
            return None

        pp = np.matrix(W)*np.matrix(P.reshape(-1,1))

        print('WESS: M*p = {}'.format(pp))
        print('WESS: max = {}'.format(np.max(pp)))

        return P, err, sactive


def prob_adjust(binprob, rates, uncert, oldindex, targets=[]):

    nbins = binprob.size

    result = solve_steady_state(rates, uncert, targets)

    if result is None:
        return binprob
    else:
        ss_estimate, err, active_bins = result

    print('WESS NNLS norm: {}'.format(err))

    #now remap new_weights onto the bin indices
    MappedActiveBins = []
    for ibin in active_bins:
        MappedActiveBins.append(oldindex[ibin])

    mapped_new_weights = np.zeros((nbins,))
    mapped_new_weights[MappedActiveBins] = ss_estimate

    # Check to make sure no bins with non-zero probability end up with zero prob after reweighting
    orig_nzi = np.nonzero(binprob)[0] # bins that originally had non-zero prob
    new_zi = np.where(mapped_new_weights == 0.0)[0] # bins that after reweight have zero prob
    ri = np.intersect1d(orig_nzi,new_zi) # bins to reset
    mapped_new_weights[ri] = 1.0E-16

    # Check to make sure no bins with zero probability end up with nonzero prob after reweighting
    orig_zi = np.where(binprob == 0.0)[0] # bins that originally had zero prob
    new_nzi = np.nonzero(mapped_new_weights)[0] # bins that after reweight have non-zero prob
    ri = np.intersect1d(orig_zi,new_nzi) # bins to reset
    mapped_new_weights[ri] = 0.0

    # Ensure that all target bins have their weight set to zero
    MappedTargetBins = []
    for ibin in targets:
        MappedTargetBins.append(oldindex[ibin])
    mapped_new_weights[MappedTargetBins] = 0.0

    mapped_new_weights /= np.sum(mapped_new_weights)

    assert not (mapped_new_weights < 0).any()

    return mapped_new_weights
