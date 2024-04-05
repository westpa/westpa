import numpy
from numpy.random import Generator, MT19937
from fasthist import histnd


def _test_double(npts=1024 * 1024, loops=3):
    from time import time

    rng = Generator(MT19937())  # RNG for the rest of this function

    mine_times = [None] * loops
    theirs_times = [None] * loops

    # though 1.0 should be sufficient, weirdness in the boundary conditions
    # for numpy.digitize appears to introduce a discrepancy, so throw something
    # greater than 1.0 in for good measure
    binbounds = [[0, 0.5, 1, 1.1] for _ in range(3)]
    weights = rng.random(size=npts)
    # weights = numpy.ones((npts,), dtype=numpy.float64)
    for n in range(loops):
        testdat = rng.random(size=(npts, 3))
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop - mstart
        theirs_times[n] = tstop - tstart
        print(mine)
        print(theirs)
        errsum = numpy.abs(mine - theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print('sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item))

    print('mine, best of {}:   {}'.format(loops, min(mine_times)))
    print('theirs, best of {}: {}'.format(loops, min(theirs_times)))


def _test_float(npts=1024 * 1024, ndim=3, loops=3):
    from time import time

    rng = Generator(MT19937())  # RNG for the rest of this function

    mine_times = [None] * loops
    theirs_times = [None] * loops

    binbounds = [[0, 0.5, 1, 1.1] for _ in range(ndim)]
    # weights = rng.random(size=npts)
    weights = numpy.ones((npts,), dtype=numpy.float64)
    for n in range(loops):
        testdat = rng.random(size=(npts, ndim), dtype=numpy.float32)
        print(testdat)
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop - mstart
        theirs_times[n] = tstop - tstart
        print(mine)
        print(theirs)
        errsum = numpy.abs(mine - theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print('sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item))

    print('mine, best of {}:   {}'.format(loops, min(mine_times)))
    print('theirs, best of {}: {}'.format(loops, min(theirs_times)))


def _test_uint(npts=1024 * 1024, ndim=3, loops=3):
    from time import time

    rng = Generator(MT19937())
    mine_times = [None] * loops
    theirs_times = [None] * loops

    binbounds = [[0, 1, 2, 3, 4] for x in range(ndim)]
    # weights = rng.random(size=npts)
    weights = numpy.ones((npts,), dtype=numpy.float64)
    print('binbounds: {}'.format(binbounds))
    print('weights')
    print(weights)
    for n in range(loops):
        testdat = rng.integers(0, 4, size=(npts, ndim), dtype=numpy.uint16)
        print('test data')
        print(testdat)
        mstart = time()
        mine = histnd(testdat, binbounds, weights=weights)
        mstop = time()
        tstart = time()
        theirs = numpy.histogramdd(testdat, binbounds, weights=weights)[0]
        tstop = time()
        mine_times[n] = mstop - mstart
        theirs_times[n] = tstop - tstart
        print(mine)
        print(theirs)
        errsum = numpy.abs(mine - theirs).sum()
        errsum_per_item = errsum / npts
        rel_err = errsum / numpy.abs(weights).sum()
        print('sum of the absolute errors: {} ({} relative, {} per entry)'.format(errsum, rel_err, errsum_per_item))

    print('mine, best of {}:   {}'.format(loops, min(mine_times)))
    print('theirs, best of {}: {}'.format(loops, min(theirs_times)))


_test_float(npts=1024 * 1024, ndim=1)
