from __future__ import division
import math
from copy import copy
import numpy

_v_floor = numpy.vectorize((lambda x: int(math.floor(x))), [numpy.int_])
_v_ceil = numpy.vectorize((lambda x: int(math.ceil(x))), [numpy.int_])

def quantile(data, q):
    N = data.shape[0]
    i_q = int(math.ceil(N*q))-1
    sdata = copy(data)
    sdata.sort()
    return sdata[i_q]

def quantiles(data, qs):
    N = data.shape[0]
    qs = numpy.array(qs)
    sdata = numpy.msort(data)
    return sdata[_v_ceil(qs*N)-1]

def quantile_range(data, q_lb, q_ub):
    N = data.shape[0]
    i_lb = int(math.floor(N * q_lb))
    i_ub = int(math.ceil(N * q_ub))
    sdata = copy(data)
    sdata.sort()
    return (sdata[i_lb], sdata[i_ub])

def ci_range(data, alpha):
    if alpha < 0:
        raise ValueError('alpha must be positive')
    elif alpha >= 0.5:
        raise ValueError('alpha must be < 0.5')
    
    halfalpha = alpha / 2
    return quantile_range(data, halfalpha, 1-halfalpha)
