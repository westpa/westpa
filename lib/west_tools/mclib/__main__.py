if __name__ == '__main__':
    from . import autocorrel_elem
    import numpy
    from scipy.signal import correlate
    
    n = 16
    x = numpy.linspace(0,n*numpy.pi,16*n+1)
    a = numpy.cos(x) + numpy.exp(-(x/2.0)**2) + numpy.exp(-(x/4.0))
    pa = numpy.zeros((10000*len(a),), numpy.float64)
    pa[:len(a)] = a
    
    print '<a> =', a.mean()
    print '<a^2> =',((a-a.mean())**2).sum()
    print 'scipy.signal.correlate:'
    acf0 = correlate(a,a)
    acf0 = acf0[-len(a):]
    acf0 /= acf0.max()
    print acf0[:len(acf0)/4]
    
#    print 'scipy.signal.correlate (-mean):'
#    acf0 = correlate(a-a.mean(),a-a.mean())
#    acf0 = acf0[-len(a):]
#    acf0 /= acf0.max()
#    print acf0[:len(acf0)/4]
        
    print 'this module:'
    acf = numpy.array([autocorrel_elem(pa,k) for k in xrange(len(a))])
    print acf[:len(acf)/4]
    

    
    
    