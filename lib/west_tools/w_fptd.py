import numpy as np
import h5py
from scipy import linalg as LA
from westpa.h5io import WESTPAH5File
from westtools import WESTDataReader, WESTTool
# Bound State: 0,2,4
# Unbound State: 1,3,5

"""
# TODO:
    implement
    go()
    process_args()
    see: core.py
    add_args()
    data manager?
    w_assign.py: WESTDataReader()
    ^see line 449

    #look at: w_assign
"""

class WFptd(WESTTool):

    def __init__(self):
        super(WFptd,self).__init__()
        self.data_reader = WESTDataReader()
        self.input_file = None
        self.output_hdf5file = None

        self.n_bins = 0
        self.init_bins = []
        self.target_bins = []
        self.cbins = []
        self.dt = 0


        self.K = None
        self.iter = 100000
        self.merged_rates = []
        self.description = 'Calculate the first passage time distribution for transitions \
            between two states.'

        # main(self) (inherited) calls make_parser_and_process() then calls go()

    def add_args(self, parser):

        self.data_reader.add_args(parser) # get HDF5 from user
        cgroup = parser.add_argument_group('Calculation options')
        cgroup.add_argument('--nbins', help='Number of bins', metavar="num_bins", dest='nbins', type=int)
        cgroup.add_argument('-i', required=True, help='Initial state bins', metavar='init_bins', type=int, nargs='+')
        cgroup.add_argument('-f', required=True, help='Final state bins', metavar='final_bins', type=int, nargs='+')
        cgroup.add_argument('--iter', required=False, help="Number of iterations to calculate the distribution over",
            metavar="iterations", type=int, nargs=1, default=1000000) # @TODO: ultimately we should be able to pull this in from the file
        cgroup.add_argument('-a', '--assign', required=False, dest='assign',
            default='assign.h5', metavar='assign_filename',
            help='Assign filename (default: \'assign.h5\')')
        ogroup = parser.add_argument_group('WEST output data options')
        # Q FOR AUDREY: what should default input/output filenames be? A FROM AUDREY: whatever you want
        ogroup.add_argument('-o', '--output', required=False, dest='output',
            default='fptd.h5', nargs=1, metavar='output_filename',
            help='Output filename (default: \'fptd.h5\')')


    def process_args(self, args):
        '''Take argparse-processed arguments associated with this component and deal
        with them appropriately (setting instance variables, etc)'''
        #args = parser.parse_args()
        self.data_reader.process_args(args) # Q: open hdf5 file for reading?
        self.input_file = h5py.File(self.data_reader.we_h5filename) # Q: is there another way to do this? A: Yes, which we should move to ultimately.
        self.output_file = WESTPAH5File(args.output, 'w', creating_program=True)
        if (args.iter):
            self.iter = args.iter
        self.dt = h5py.File(args.assign)['assignments'].shape[2]
        self.set_bin_info(args)
        self.set_K(args)

    def set_bin_info(self, args):
        self.n_bins = args.nbins
        self.init_bins = args.i
        self.target_bins = args.f
        # cbins: bins in neither initial nor target state
        self.cbins = list(x for x in range(0, self.n_bins) if x not in self.target_bins)
        self.merged_rates = np.empty(self.n_bins - len(self.target_bins)) # size: num bins which are NOT Target bins

    def set_K(self, args):
        f = self.input_file # convenience variable
        # Empty shell for sparse trans matrix with weights.
        self.K = np.zeros((self.n_bins,self.n_bins))
        #trans_m = self.K.copy()
        # Fill matrix with flux vals from h5 file.
        for iter in f['iterations'].keys():
            grp = f['iterations'][iter]
            rows = grp['rows']
            cols = grp['cols']
            flux = grp['flux']
            self.K[rows, cols] += flux
        # Normalize flux values
        #self.K = trans_m
        for row in range(self.n_bins):
            #if self.K[row, :].sum() == 0:
            #    dell.append(row) # AP: dell?
            #else:
            self.K[row, :] /= self.K[row, :].sum()
        # NaN values are set to zero, infinity values are set to large finite numbers.
        self.K = np.nan_to_num(self.K)

    def go(self):
        '''Perform the analysis associated with this tool.'''
        # eigenvalues, eigenvectors of tranpose of K
        # eigenvalues: solutions to characteristic polynomial for K
        # eigenvectors: solutions v to eqn. K<dot>v = lambda<dot>v
        eigvals, eigvecs = LA.eig(self.K.T)
        unity = (np.abs(np.real(eigvals) - 1)).argmin()
        # print("eigenvalues", eigvals)
        # print("eigenvectors", eigvecs)
        # print("unity", unity)
        eq_pop = np.abs(np.real(eigvecs)[unity])
        eq_pop /= eq_pop.sum()

        # probability of starting in init bin A.
        # distr_prob = np.random.rand(len(INIT_BINS))
        paths = []
        # t_bins: all bins which are not target bins.
        # t_bins = list(x for x in range(0, NBINS) if x not in TARGET_BINS)
        # lower_bound = mfpt - error
        # lower_bound = 121.8
        # lower_bound = 116 # AP: where did the value for lower_bound come from?  NOT NECESSARY

        # delete (?) (lines 141-144 in new_paths_annotated):
        """
        p_dist = np.zeros(NBINS)
        p_dist[INIT_BINS] = 1.0/float(len(INIT_BINS))
        pp_dist = np.zeros(NBINS)
        p_dist = eq_pop
        p_dist = p_dist
        p_dist = np.zeros((NBINS,NBINS))
        """
        p_dist = self.K.copy()

        histogram = []
        # is multiplying by the eq_pop sum == mutiplying by probability of starting in given init bin?
        # do we need to add this line below?
        histogram.append(p_dist[self.init_bins, self.target_bins]*eq_pop[self.cbins].sum())
        # bc the first iteration below would start adding to histogram at equiv of n=2 steps
        #ITER = 1600000 # AP: should this be self.iter?  Cool, seems like 'yes'.
        for i in range(self.iter):
            # np.diag(...) : gives 1D arr of diag vals of matrix
            # np.diag(np.diag(...)) : gives the actual diagonal of matrix
            np_dist = np.dot(self.K, p_dist - np.diag(np.diag(p_dist)))
            histogram.append(np_dist[self.init_bins, self.target_bins]*eq_pop[self.cbins].sum())
            p_dist = np_dist

        #dt = 101 # AP: did this come from a specific simulation?

        print(np.nan_to_num(histogram).shape, len(range(0, self.iter+1)))
        # MFPT?
        print(self.iter, (np.average(range(0,self.iter+1), weights=np.nan_to_num(histogram)[:,0])/self.dt))
        print(eq_pop)
        print(eq_pop[self.cbins].sum())
        # for graph: would sum each thing in histogram?
        # Save results
        # "histogram"
        # in paper: prob density vs FPT(A->B)/ns
        # FPT - x axis, prob density - y axis
        self.output_file.create_dataset('probability_densities', dtype=np.float64, shape=())
        self.output_file.create_dataset('initial_states')
        self.output_file.create_dataset('final_states')

if __name__ == '__main__':
    WFptd().main()
