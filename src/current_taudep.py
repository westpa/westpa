from __future__ import division
import sys
from optparse import OptionParser
import numpy
import matplotlib
from math import sqrt, floor, ceil

matplotlib.rcParams['backend'] = 'PDF'
from matplotlib import pyplot


def mc_bootstrap(data, n_sets):
    synth = numpy.empty((n_sets,) + data.shape, data.dtype)
    N = data.shape[0]
    for i in xrange(0, n_sets):
        indices = numpy.random.randint(N,size=(N,))
        synth[i] = data[indices]
    return synth

def mc_bootstrap_ci(estimator, data, alpha, n_sets=1000):
    alpha = min(alpha, abs(1-alpha))
    half_alpha = alpha / 2.0

    fhat = estimator(data)
    f_synth = numpy.empty((n_sets,), data.dtype)
    synth = mc_bootstrap(data, n_sets)
    for i in xrange(0, n_sets):
        f_synth[i] = estimator(synth[i])
    f_synth.sort()
    lb = f_synth[int(floor(n_sets*half_alpha))]
    ub = f_synth[int(ceil(n_sets*(1-half_alpha)))]

    return (fhat, lb, ub)

def blockavg(data,
             output_filename = None, 
             plot_filename = None,
             title = 'Block Average Error Analysis',
             xlabel = 'Number of blocking transformations',
             ylabel = r'$\hat\sigma_{\left<x\right>}\ \pm \sigma_{\hat\sigma}$'):
    
    stds = []
    working_data = data.copy()
    while len(working_data) > 2:
        npoints = len(working_data)
        working_mean = working_data.mean()

        varest = ((working_data-working_mean)**2).sum() / (npoints*(npoints-1))
        stdest = sqrt(varest)
        stdstd = stdest/sqrt(2*(npoints-1))

        stds.append((stdest,stdstd))

        if len(working_data) % 2 == 1:
            working_data[-2] = (working_data[-2] + working_data[-1]) / 2
            working_data = working_data[:-1]
        working_data = (working_data[::2] + working_data[1::2]) / 2
    
    stds = numpy.array(stds)
    if output_filename:
        numpy.savetxt(output_filename, stds)

    if plot_filename:
        pyplot.errorbar(range(0,len(stds)), stds[:,0], yerr=stds[:,1],
                        zorder=10)
        pyplot.axis([0,len(stds),None,None])
        if title:
            pyplot.title(title)
        pyplot.xlabel(xlabel)
        pyplot.ylabel(ylabel)
        pyplot.savefig(plot_filename)
        pyplot.clf()

    for itrans in xrange(len(stds)-1, 1, -1):
        lb = stds[itrans,0] - stds[itrans,1]
        ub = stds[itrans,0] + stds[itrans,1]
        if not (lb <= stds[itrans-1,0] < ub):
            break
    else:
        print "block average did not converge"
        itrans = -1


    return (itrans, stds[itrans,0], stds[itrans,1], stds)

parser = OptionParser(usage='%prog [OPTIONS] CURRENT_FILE')
parser.add_option('-b', '--start', '--begin', dest='start', type='int',
                  help='begin analysis at START (default: 0)',
                  default=0)
parser.add_option('-e', '--stop', '--end', dest='stop', type='int',
                  help='end analysis at STOP (inclusive; default: all data)')
parser.add_option('-s', '--step', dest='step', type='int',
                  help='output analysis for every STEP current values '
                      +'(default: 100)',
                  default=100)
parser.add_option('--skip-zero', dest='skip_zero', type='int',
                  help='skip SKIP_ZERO times the number of zero current '
                      +'values (default: 2)',
                  default=2)
parser.add_option('--skip-points', dest='skip_points', type='int',
                  help='skip SKIP_POINTS (default: 0)',
                  default=0)
parser.add_option('--skip-fraction', dest='skip_fraction', type='float',
                  help='skip SKIP_FRACTION of the data (default 0.0)',
                  default=0.0)
parser.add_option('--blockavg-set-output', dest='blockavg_set_output',
                  help='record each set of block averaging results in '
                      +'BLOCKAVG_SET_OUTPUT; may contain a %d substitution '
                      +'for N_tau (default: current_blockavg_%d.txt)',
                  default='current_blockavg_%d.txt')
parser.add_option('--blockavg-set-plot', dest='blockavg_set_plot',
                  help='plot each set of block averaging results to '
                      +'BLOCKAVG_SET_PLOT; '
                      +'may contain a %d substitution for N_tau '
                      +'(default: current_blockavg_%d.pdf)',
                  default='current_blockavg_%d.pdf')
parser.add_option('--mc-alpha', dest='mc_alpha', type='float',
                  help='significance level for Monte Carlo bootstrap',
                  default=0.05)
parser.add_option('--mc-nsynth', dest='mc_nsynth', type='int',
                  help='number of synthetic data sets for MC bootstrap '
                      +'(default: 2000)',
                  default=2000)
parser.add_option('--blockavg-output', dest='blockavg_output',
                  help='store block average analysis results in '
                      +'BLOCKAVG_OUTPUT (default: current_taudep_ba.txt)',
                  default='current_taudep_ba.txt')
parser.add_option('--blockavg-plot', dest='blockavg_plot',
                  help='plot block average analysis results to '
                      +'BLOCKAVG_PLOT (default: current_taudep_ba.pdf)',
                  default='current_taudep_ba.pdf')
parser.add_option('--mc-output', dest='mc_output',
                  help='store MC bootstrap analysis results in MC_OUTPUT '
                      +'(default: current_taudep_mc.txt)',
                  default='current_taudep_mc.txt')
parser.add_option('--mc-plot', dest='mc_plot',
                  help='plot MC bootstrap analysis results to MC_PLOT '
                      +'(default: current_taudep_mc.pdf',
                  default='current_taudep_mc.pdf')

(opts, args) = parser.parse_args()

if not args:
    sys.stderr.write('a file containing current values is required\n')
    parser.print_help(sys.stderr)
    sys.exit(1)
else:
    print "loading current values from '%s'" % args[0]

current = numpy.loadtxt(args[0])
numzero = len(current[current[:,2] == 0.0])
numtrunc = max([int(opts.skip_zero) * numzero,
                int(opts.skip_points),
                int(opts.skip_fraction * current.shape[0])])

print "omitting first %d points" % numtrunc

td_current_ba = []
td_current_mc = []
start = opts.start
if opts.stop:
    stop = opts.stop
else:
    stop = len(current)+1
step = opts.step
for n_tau in xrange(start, stop, step):
    if n_tau <= numtrunc: continue

    if opts.blockavg_set_output:
        try:
            bavg_output = opts.blockavg_set_output % n_tau
        except TypeError:
            bavg_output = '%s_%d' % (opts.blockavg_set_output, n_tau)
    else:
        bavg_output = None

    if opts.blockavg_set_plot:
        try:
            bavg_plot = opts.blockavg_set_plot % n_tau
        except TypeError:
            bavg_plot = '%s_%d' % (opts.blockavg_set_plot, n_tau)
    else:
        bavg_plot = None

    print "\nn_tau = %d" % n_tau
    avg_current = current[numtrunc:n_tau,2].mean()
    bavg_data = blockavg(current[numtrunc:n_tau, 2],
                         bavg_output, bavg_plot,
                         title = r'Block Average Error Analysis ($N_\tau=%d$)'%n_tau)
    (ntrans, sem, sem_error, stds) = bavg_data
    print "  block average results:"
    print "    correlation time:    %d tau" % ntrans
    print "    block average S.E.M: %g" % sem
    print "    rate +/- 2 S.E.M.:   %g +/- %g" % (avg_current, 2*sem)

    td_current_ba.append((n_tau, avg_current, 2*sem))

    (mc_avg, mc_lb, mc_ub) = mc_bootstrap_ci(numpy.mean, 
                                             current[numtrunc:n_tau, 2], 
                                             opts.mc_alpha, 
                                             opts.mc_nsynth)
    print "  MC bootstrap results:"
    print "    95%% CI lower bound:  %g" % mc_lb
    print "    average:             %g" % mc_avg
    print "    95%% CI upper bound:  %g" % mc_ub
    print "    95%% CI width:        %g" % (mc_ub - mc_lb)

    td_current_mc.append((n_tau, mc_avg, mc_lb, mc_ub))


if opts.blockavg_output:
    numpy.savetxt(opts.blockavg_output, numpy.array(td_current_ba))

if opts.mc_output:
    numpy.savetxt(opts.mc_output, numpy.array(td_current_mc))

if opts.blockavg_plot:
    pyplot.clf()
    pyplot.errorbar([tdc[0] for tdc in td_current_ba],
                    [tdc[1] for tdc in td_current_ba],
                    [tdc[2]*2 for tdc in td_current_ba])
    pyplot.ylabel(r'$\left<k\right>\pm2\sigma(\left<k\right>)$')
    pyplot.xlabel(r'$N_\tau$')
    pyplot.title('Current Analysis (Block Averaging)')
    pyplot.savefig(opts.blockavg_plot)

if opts.mc_plot:
    pyplot.clf()
    td_current_bounds = numpy.empty((len(td_current_mc), 2), numpy.float64)
    for i in xrange(0, len(td_current_mc)):
        td_current_bounds[i,0] = td_current_mc[i][2]
        td_current_bounds[i,1] = td_current_mc[i][3]
    pyplot.errorbar([tdc[0] for tdc in td_current_mc],
                    [tdc[1] for tdc in td_current_mc],
                    td_current_bounds.T)
    pyplot.ylabel(r'$\left<k\right>$ and 95% CI')
    pyplot.xlabel(r'$N_\tau$')
    pyplot.title('Current Analysis (Monte Carlo Bootstrap)')
    pyplot.savefig(opts.mc_plot)

