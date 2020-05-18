
import numpy
from mclib import mcbs_ci, mcbs_ci_correl, mcbs_correltime

n_sets = 1000

#N = 10000
#correl_len = 16
#ds_uncorrel = numpy.array(numpy.random.random_integers(0,1, size=(N,)), dtype=numpy.float64)
#ds_correl = numpy.empty((N,), dtype=numpy.float64)
#for i in xrange(0,N,correl_len):
#    ds_correl[i:i+correl_len] = numpy.random.randint(0,2)

ds_uncorrel = ds_correl = numpy.loadtxt('/home/mzwier/flexible-p53-MDM2-1d-rmsd/avgflux_1_blocked.txt', usecols=[2])

    
print('mean:')
print('  uncorrelated:    {}'.format(ds_uncorrel.mean()))
print('  correlated:      {}'.format(ds_correl.mean()))

print('correlation time:')
k_uncorrel = mcbs_correltime(ds_uncorrel,0.05,n_sets)
k_correl = mcbs_correltime(ds_correl,0.05,n_sets)
print('  uncorrelated:    {}'.format(k_uncorrel))
print('  correlated:      {}'.format(k_correl))


print('naive MCBS:')
uncorrel_mean, uncorrel_lb, uncorrel_ub = mcbs_ci(ds_uncorrel, numpy.mean, 0.05, n_sets=n_sets, sort=numpy.msort)
correl_mean, correl_lb, correl_ub = mcbs_ci(ds_correl, numpy.mean, 0.05, n_sets=n_sets, sort=numpy.msort)
n_uncorrel_width = uncorrel_ub - uncorrel_lb
n_correl_width   = correl_ub - correl_lb
print('  uncorrelated:    {} ({},{})'.format(uncorrel_mean, uncorrel_lb,uncorrel_ub))
print('  correlated:      {} ({},{})'.format(correl_mean, correl_lb,correl_ub))
print('  width ratio c/u: {}'.format((n_correl_width/n_uncorrel_width)))

print('blocked MCBS:')
#uncorrel_mean, uncorrel_lb, uncorrel_ub = mcbs_ci(ds_uncorrel[::k_uncorrel+1], numpy.mean, 0.05, n_sets=1000, sort=numpy.msort)
#correl_mean, correl_lb, correl_ub = mcbs_ci(ds_correl[::k_correl+1], numpy.mean, 0.05, n_sets=1000, sort=numpy.msort)
uncorrel_mean, uncorrel_lb, uncorrel_ub, k_ = mcbs_ci_correl(ds_uncorrel, numpy.mean, 0.05, n_sets=n_sets, subsample=numpy.mean)
correl_mean, correl_lb, correl_ub, k_ = mcbs_ci_correl(ds_correl, numpy.mean, 0.05, n_sets=n_sets, subsample=numpy.mean)
b_uncorrel_width = uncorrel_ub - uncorrel_lb
b_correl_width   = correl_ub - correl_lb
print('  uncorrelated:    {} ({},{})'.format(uncorrel_mean, uncorrel_lb,uncorrel_ub))
print('  correlated:      {} ({},{})'.format(correl_mean, correl_lb,correl_ub))
print('  width ratio c/u: {}'.format((b_correl_width/b_uncorrel_width)))

print('width ratio blocked/naive:')
print('  uncorrelated:    {}'.format(b_uncorrel_width/n_uncorrel_width))
print('  correlated:      {}'.format(b_correl_width/n_correl_width))

