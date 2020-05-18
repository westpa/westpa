import argparse
import h5py
from matplotlib import pyplot

parser = argparse.ArgumentParser(description='''Plot flux evolution obtained with w_fluxanl.''')
parser.add_argument('-s', '--state', default=0,
                    help='''Create a plot for target state index STATE (default: %(default)s)''')
parser.add_argument('-i', '--input', default='fluxanl.h5',
                    help='''Take data from INPUT (default: %(default)s)''')
parser.add_argument('-o', '--output', default='fluxevol.pdf',
                    help='''Write plot of flux evolution to OUTPUT (default: %(default)s)''')
args = parser.parse_args()

h5file = h5py.File(args.input, 'r')
evol = h5file['flux_evol'][...]
evol_iters = h5file['flux_evol_iterations'][...]

pyplot.figure()
iters = evol_iters['iter_stop']-1
pyplot.plot(iters, evol['mean'], color='black')
pyplot.plot(iters, evol['ci_lb'], color='gray')
pyplot.plot(iters, evol['ci_ub'], color='gray')
pyplot.title(r'Flux evolution for state {} at $\alpha=${} confidence interval'
             .format(args.state, h5file['flux_evol'].attrs['mcbs_alpha']))
pyplot.xlabel('Iteration')
pyplot.ylabel(r'Flux ($\tau^{-1}$)')
pyplot.savefig(args.output)
