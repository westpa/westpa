import numpy as np
import h5py
from scipy import linalg as LA
from westpa.h5io import WESTPAH5File
from westtools import WESTDataReader
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

"""

class WFptd(WESTTool):

    def __init__(self):
        super(WFptd,self).__init__()
        self.data_reader = WESTDataReader()
        self.n_bins = 0
        self.init_bins = []
        self.target_bins = []
        self.cbins = []
        self.K = None
        self.merged_rates = []
        self.description = 'Calculate the first passage time distribution for transitions \
            between two states.'
        make_parser()
        add_args(parser)

    # @ self, parser
    def add_args(self, parser):
        self.data_reader.add_args(parser) # get HDF5 from user
        parser.add_argument('nbins', help='Number of bins', metavar="num_bins", type=int)
        cgroup = parser.add_argument_group('Calculation options')
        cgroup.add_argument('-i', required=True, help='Initial state bins', metavar='init_bins', type=int, nargs='+')
        cgroup.add_argument('-f', required=True, help='Final state bins', metavar='final_bins', type=int, nargs='+')

        # check if sum of init/final state bin is <= total # of bins?
        #
        # output option?
    # @ self, args
    def process_args(self, args):
        '''Take argparse-processed arguments associated with this component and deal
        with them appropriately (setting instance variables, etc)'''
        self.data_reader.open(mode='r')
        args = parser.parse_args()
        #with self.data_reader: # opens HDF5 file?
            #if args.config_from_file == False:
                #continue
                    #self.binning.set_we_h5file_info(self.n_iter,self.data_reader)
                    #self.binning.process_args(args)

            #self.output_filename = args.output
        #set_bin_info()
        # set up K? (transition matrix)

    def set_bin_info(n_bins, init_bins, target_bins):
        self.n_bins = n_bins
        self.init_bins = init_bins
        self.target_bins = target_bins
        # cbins: bins in neither initial nor target state
        self.cbins = list(x for x in range(0, self.n_bins) if x not in self.target_bins)
        self.merged_rates = np.empty(self.n_bins - len(self.target_bins)) # size: num bins which are NOT Target bins

    def go(self):
        '''Perform the analysis associated with this tool.'''
        # 1. Get cbins: bins not in
        # 2.
        # 3. get transition matrix K
        # 4. open h5 file
        # 5. trans_m: empty shell for sparse trans matrix w/weights
        # 6. fill with flux info
        dell = []
        for row in range(self.n_bins):
            if trans_m[row, :].sum() == 0:
                dell.append(row)
            else:
                trans_m[row, :] /= trans_m[row, :].sum()

        # NaN values are set to zero, infinity values are set to large finite numbers.
        K = np.nan_to_num(trans_m)
        print(K)

        # eigenvalues, eigenvectors of transpose of K
        eigvals, eigvecs = LA.eig(K.T)
        unity = (np.abs(np.real(eigvals) - 1)).argmin()  # index of Smallest value from (abs of (real part of eigenvalues - 1))
        print(eigvals, eigvecs)
        eq_pop = np.abs(np.real(eigvecs)[unity])  # abs val of real part of eigenvec at index of unity
        eq_pop /= eq_pop.sum()  # divide by self? what is this
        # TODO: change this to actual distr_prob
        # probability of starting in init bin A.
        distr_prob = np.random.rand(len(self.init_bins))
        paths = []

        # lower_bound = mfpt - error
        # lower_bound = 121.8
        lower_bound = 116


        # Hmmm.  What about...
        p_dist = np.zeros(self.n_bins)
        p_dist[self.init_bins] = 1.0/float(len(self.init_bins))
        pp_dist = np.zeros(self.n_bins)
        p_dist = eq_pop

        p_dist = p_dist
        p_dist = np.zeros((self.n_bins, self.n_bins))
        p_dist = K.copy()

        histogram = []

        ITER = 1600000
        # ITER = 0
        # 100 iterations?
        for i in range(ITER):
            np_dist = np.dot(K,
                             p_dist - np.diag(np.diag(p_dist)))
            histogram.append(np_dist[self.init_bins, self.target_bins]*eq_pop[self.cbins].sum())
            p_dist = np_dist

        dt = 101

        print(np.nan_to_num(histogram).shape, len(range(1, ITER+1)))
        print(ITER, (np.average(range(1,ITER+1), weights=np.nan_to_num(histogram)[:,0])/dt))
        print(eq_pop)
        print(eq_pop[self.cbins].sum())

if __name__ == '__main__':
    WFptd().main()
