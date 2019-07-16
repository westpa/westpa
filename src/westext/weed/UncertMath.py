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
import itertools

TOL=1.0E-6

class UncertContainer(object):
    """ Container to hold uncertainty measurements. Data is convert to np masked arrays
        to avoid possible numerical problems
    """
    def __init__(self,vals,vals_dmin,vals_dmax,mask=ma.nomask):
        super(UncertContainer, self).__init__()
        
        # If input data already masked arrays extract unmasked data
        if ma.isMaskedArray(vals): 
            vals = vals.data
        if ma.isMaskedArray(vals_dmin):
            vals_dmin = vals_dmin.data
        if ma.isMaskedArray(vals_dmax):
            vals_dmax = vals_dmax.data
        
        # Adjust negative values
        ineg = np.where(vals_dmin <= 0.0)
        vals_dmin[ineg] = TOL*vals[ineg]

        # Calculate weight based on fractional uncertainty 
        diff = vals_dmax - vals_dmin
        diff_m = ma.masked_where(vals_dmax == vals_dmin,diff)        

        self.vals = ma.masked_where(vals == 0.0,vals)

        self.wt = (self.vals/diff_m)**2
        self.uncert = diff_m/self.vals

        self.wt.fill_value = np.inf
        self.uncert.fill_vaule = np.inf

        assert np.all(self.wt.mask == self.uncert.mask)
        
        # Mask data if uncertainty is not finite or if any of the inputs were
        # already masked

        mm = ma.mask_or(self.wt.mask,mask)
        
        self.vals.mask = mm
        self.wt.mask = mm
        self.uncert.mask = mm
        self.dmin = ma.array(vals_dmin,mask=mm,fill_value=np.inf)
        self.dmax = ma.array(vals_dmax,mask=mm,fill_value=np.inf)

        self.mask = ma.getmaskarray(self.vals)

    def __getitem__(self,indx):
        vals = self.vals[indx]
        dmin = self.dmin[indx]
        dmax = self.dmax[indx]
        
        if isinstance(vals, ma.core.MaskedConstant):
            dum = np.zeros((1,))
            return UncertContainer(dum.copy(),dum.copy(),dum.copy())
        elif isinstance(vals,(float,int,np.float,np.int)):
            return UncertContainer(np.array([vals,]),np.array([dmin,]),np.array([dmax,]))
        elif isinstance(vals,np.ndarray):
            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        else:
            raise TypeError

    def __setitem__(self,indx,value):
        if isinstance(value,UncertContainer):
            self.vals[indx] = value.vals
            self.dmin[indx] = value.dmin
            self.dmax[indx] = value.dmax
            self.wt[indx] = value.wt
            self.uncert[indx] = value.uncert
        else:
            raise TypeError('Can only set values with an UncertContainer object')


    def __repr__(self):
        return 'shape={} vals={} dmin={} dmax={} vals.mask={}'.format(self.vals.shape,self.vals,self.dmin,self.dmax,self.vals.mask)

    def __add__(self,value):
        if isinstance(value,UncertContainer):
            vals = self.vals + value.vals
            dmin = self.dmin + value.dmin
            dmax = self.dmax + value.dmax

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        elif isinstance(value,(float,int,np.float,np.int)):
            vals = self.vals + value
            dmin = self.dmin + value
            dmax = self.dmax + value

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        else:
            raise TypeError('Attempt to add value of unsupported type')

    def __sub__(self,value):
        if isinstance(value,UncertContainer):
            vals = self.vals - value.vals
            dmin = self.dmin - value.dmin
            dmax = self.dmax - value.dmax

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        else:
            raise TypeError ('Attempted to subtract by value of unsupported type')

    def __mul__(self,value):
        if isinstance(value,UncertContainer):
            vals = self.vals * value.vals
            dmin = self.dmin * value.dmin
            dmax = self.dmax * value.dmax

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        
        elif isinstance(value,(float,int,np.float,np.int)):
            vals = self.vals * value
            dmin = self.dmin * value
            dmax = self.dmax * value

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        else:
            raise TypeError('Attempted to multiply by value of unsupported type')

    def __div__(self,value):
        if isinstance(value,UncertContainer):
            vals = self.vals / value.vals
            dmin = self.dmin / value.dmax
            dmax = self.dmax / value.dmin

            return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        else:
            raise TypeError('Attempted to divide by unsupported type')

    def transpose(self):
        vals = self.vals.T
        dmin = self.dmin.T
        dmax = self.dmax.T

        return UncertContainer(vals,dmin,dmax,mask=vals.mask)
        
    def recip(self):
        vals = 1.0 / self.vals
        dmin = 1.0 / self.dmax
        dmax = 1.0 / self.dmin

        return UncertContainer(vals,dmin,dmax,mask=vals.mask)

    def update_mask(self):
        self.vals.mask = self.mask
        self.dmin.mask = self.mask
        self.dmax.mask = self.mask
        self.wt.mask = self.mask
        self.uncert.mask = self.mask

    def concatenate(self,value,axis=0):
        """ Concatentate UncertContainer value to self.
            Assumes that if dimensions of self and value do not match, to 
            add a np.newaxis along axis of value
        """

        if isinstance(value,UncertContainer):
            if value.vals.ndim == self.vals.ndim:
                vals = value.vals
                dmin = value.dmin
                dmax = value.dmax
                wt = value.wt
                uncert = value.uncert
                mask = value.mask
            elif (value.vals.ndim + 1) == self.vals.ndim:
                vals =  ma.expand_dims(value.vals,axis)
                dmin =  ma.expand_dims(value.dmin,axis)
                dmax =  ma.expand_dims(value.dmax,axis)
                wt =  ma.expand_dims(value.wt,axis)
                uncert =  ma.expand_dims(value.uncert,axis)
                mask =  np.expand_dims(value.mask,axis)
            else:
                raise ValueError('Could not propery match dimensionality')
                
            self.vals = ma.concatenate((self.vals,vals),axis=axis)
            self.dmin = ma.concatenate((self.dmin,dmin),axis=axis)
            self.dmax = ma.concatenate((self.dmax,dmax),axis=axis)
            self.wt = ma.concatenate((self.wt,wt),axis=axis)
            self.uncert = ma.concatenate((self.uncert,uncert),axis=axis)
            
            self.mask = np.concatenate((self.mask,mask),axis=axis)
        else:
            raise ValueError('Can only concatenate with an UncertContainer object')

    def weighted_average(self,axis=0,expaxis=None):
        """ Calculate weighted average of data along axis
            after optionally inserting a new dimension into the
            shape array at position expaxis
        """

        if expaxis is not None:
            vals = ma.expand_dims(self.vals,expaxis)
            dmin = ma.expand_dims(self.dmin,expaxis)
            dmax = ma.expand_dims(self.dmax,expaxis)
            wt = ma.expand_dims(self.wt,expaxis)
        else:
            vals = self.vals
            wt = self.wt
            dmin = self.dmin
            dmax = self.dmax
        
        # Get average value
        avg,norm = ma.average(vals,axis=axis,weights=wt,returned=True)
        avg_ex = ma.expand_dims(avg,0)

        # Calculate weighted uncertainty
        wtmax = ma.max(wt,axis=axis)
        neff = norm/wtmax       # Effective number of samples based on uncertainties

        # Seeking max deviation from the average; if above avg use max, if below use min
        term = np.empty_like(vals)
        
        indices = np.where(vals > avg_ex)
        i0 = indices[0]
        irest = indices[1:]
        ii = tuple(x for x in itertools.chain([i0],irest))
        jj = tuple(x for x in itertools.chain([np.zeros_like(i0)],irest))
        term[ii] = (dmax[ii] - avg_ex[jj])**2
        
        indices = np.where(vals <= avg_ex)
        i0 = indices[0]
        irest = indices[1:]
        ii = tuple(x for x in itertools.chain([i0],irest))
        jj = tuple(x for x in itertools.chain([np.zeros_like(i0)],irest))
        term[ii] = (avg_ex[jj] - dmin[ii])**2
        
        dsum = ma.sum(term*wt,axis=0)     # Sum for weighted average of deviations

        dev = 0.5*np.sqrt(dsum/(norm*neff))
        
        if isinstance(avg,(float,np.float)):
            avg = avg_ex

        tmp_min = avg - dev
        ii = np.where(tmp_min < 0)
        tmp_min[ii] = TOL*avg[ii]
        
        return UncertContainer(avg,tmp_min,avg+dev)


