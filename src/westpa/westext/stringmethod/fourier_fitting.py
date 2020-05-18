
import numpy as np
import scipy
import scipy.optimize


class FourierFit(object):
    def __init__(self, P=2, ndims=2, maxiters=100, tol=1.0E-6):
        super(FourierFit, self).__init__()
        
        self.P = P
        self.maxiters = maxiters
        self.ndims = ndims
        self.tol = tol
        self.pp = []
        self.t0 = None
        self.w0 = None
        
    def calc_string(self,w,t,x_meas):
        tlen = len(t)
        t = np.linspace(0.0,1.0,tlen)
        x_est = x_meas[0,:] + (x_meas[-1,:] - x_meas[0,:])*t[:,np.newaxis]
        for i in range(self.ndims):
            for j in range(self.P):
                x_est[:,i] += w[i,j]*np.sin((j+1)*np.pi*t)
        return x_est

    def _optimize_dist(self,tk,x_meas,w,k):
        x_target = x_meas[k,:]
        x_est = x_meas[0,:] + (x_meas[-1,:] - x_meas[0,:])*tk
        for i in range(self.ndims):
            for j in range(self.P):
                x_est[i] += w[i,j]*np.sin((j+1)*np.pi*tk)

        err = x_target - x_est
        return err

    def _optimize_w(self,w,x_meas,t,k,weight):
        x_target = x_meas[:,k]
        x_est = x_meas[0,k] + (x_meas[-1,k] - x_meas[0,k])*t

        for j in range(self.P):
            x_est += w[j]*np.sin((j+1)*np.pi*t)

        err = weight*(x_target - x_est)
        return err
        
    def optimize(self,data,weight,w0,t0):
        ncenters = data.shape[0]
        self.w0 = w0
        self.t0 = t0
        if weight is None:
            weight = np.ones_like(t0)

        for iiter in range(self.maxiters):
            self.pp.append(self.calc_string(self.w0,self.t0,data))
            if iiter > 0:
                err = np.sum((self.pp[-1] - self.pp[-2])**2)/ncenters
                print('{} -- {}'.format(iiter,err))
                if err < self.tol:
                    break
            else:
                print(iiter)
            # Optimize tk
            for ci in range(ncenters):
                self.t0[ci] = scipy.optimize.leastsq(self._optimize_dist, self.t0[ci], args=(data,self.w0,ci))[0]
            
            # Optimize wij
            for k in range(self.ndims):
                self.w0[k,:] = scipy.optimize.leastsq(self._optimize_w,self.w0[k,:],args=(data,self.t0,k,weight))[0]
