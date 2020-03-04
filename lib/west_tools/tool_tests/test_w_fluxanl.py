
import nose
import nose.tools
from nose.tools import with_setup
from nose.plugins.skip import SkipTest

import tempfile, math
from .common import *
from w_fluxanl import WFluxanlTool
import h5py


def get_flux_group(f):
    '''f is an h5py file type object'''

    return f['target_flux']['index']

class Test_W_Fluxanl_Args(CommonToolTest):
    '''
    This class tests w_fluxanl functions with different combinations of command line arguments'''

    arg_combos = {'alpha':[None, 0.01, 0.20], 'autocorrel_alpha':[None, 0.01], 'nsets':[None, 500],
                  '--evol':[None,True], 'evol_step':[None]}

    mean_flux = None

    test_name = 'W_FLUXANL'

    def test_args(self):
        '''Testing arg combos: w_fluxanl runs as expected'''
        
        for args_dict in cycle_args(self.arg_combos):
            self.outfile = self.mktemp(prefix='fluxanl')
            self.w = WFluxanlTool()
            args = make_args(self.outfile, **args_dict)
            self.w.make_parser_and_process(args=args)

            #test_outcome, err = self.check_args_processed(**args_dict)
            test_outcome, err = self.check_runs_with_args(**args_dict)

            yield self.check_args, test_outcome, err, args[:-1]

            del self.w, self.outfile

    def check_output(self, **kwargs):
        '''Check that output file for w_fluxanl is set up correctly'''

        with h5py.File(self.outfile) as f:
            assert 'target_flux' in list(f.keys()), "'target flux' group not in output file"
            target_group = f['target_flux']
            assert 'index' in list(target_group.keys()), "'index' group not in output file"
            index_group = target_group['index']

            mean_flux = index_group['mean_flux'][0]

            if self.mean_flux is not None:
                assert mean_flux == self.mean_flux, "Actual mean flux ({}) is different than expected ({})".format(mean_flux,self.mean_flux)
            else:
                self.mean_flux = mean_flux

            
            expected_alpha = kwargs['alpha'] or 0.05
            actual_alpha = index_group.attrs['mcbs_alpha']

            assert expected_alpha == actual_alpha, "Actual alpha ({}) does not match expected ({})".format(actual_alpha, expected_alpha)

            expected_acalpha = kwargs['autocorrel_alpha'] or expected_alpha
            actual_acalpha = index_group.attrs['mcbs_autocorrel_alpha']

            assert expected_acalpha == actual_acalpha, "Actual autocorrel alpha ({}) does not match expected ({})".format(actual_alpha, expected_alpha)

            expected_nsets = kwargs['nsets'] or 10**(math.ceil(-math.log10(expected_alpha)) + 1)
            actual_nsets = index_group.attrs['mcbs_n_sets']

            assert expected_nsets == actual_nsets, "Actual n_sets ({}) does not match expected ({})".format(actual_nsets,expected_nsets)

    '''Unused method - rather than check instance variables, check output is as expected
    def check_args_processed(self,alpha=None,autocorrel_alpha=None,nsets=None,evol=False,evol_step=None):
        try:
            wtool = self.w

            #Check that dependencies, namely iter_range and data_reader objects are initialized correctly
            iter_range = wtool.iter_range
            data_reader = wtool.data_reader

            assert data_reader and iter_range, 'Data reader not initialized'

            iter_start = iter_range.iter_start
            iter_stop = iter_range.iter_stop or data_reader.current_iteration
            n_iter = iter_stop - iter_start

            assert iter_start == 1, "Iter_reader does not start at first iteration by default"
            assert n_iter == 30, "Total number of iterations is {} (expected 20)".format(n_iter)

            ##Check that fluxanl specific instance variables are set correctly
            assert wtool.alpha == (alpha or 0.05), 'Alpha set to {} (expected {})'.format(wtool.alpha, alpha or 0.05)
            assert wtool.autocorrel_alpha == (autocorrel_alpha or wtool.alpha), "acAlpha set to {} (expected {})".format(wtool.autocorrel_alpha, autocorrel_alpha or wtool.alpha)

            nsets = nsets or 10**(math.ceil(-math.log10(wtool.alpha)) + 1)
            assert wtool.n_sets == nsets, "N_sets set to {} (expected {})".format(wtool.n_sets, nsets)

            if evol:
                assert wtool.do_evol, "Evol should have been set to True"
                assert wtool.evol_step == evol_step or 1, "Evol step set to {}, expected {}".format(wtool.evol_step, evol_step or 1)

        except AssertionError as e:
            return (0, e)

        else:
            return (1, None)
    '''

@SkipTest
class Test_W_Fluxanl_System(CommonToolTest):
    '''System level tests for w_fluxanl'''

    def setUp(self):
        self.w = WFluxanlTool()
        self.outfile = self.mktemp(prefix='fluxanl')

    def tearDown(self):
        del self.w

    def test_alpha(self):
        '''Confidence interval size decreases as alpha increases'''

        default_alpha_tester = self.w #alpha == 0.05
        low_alpha_tester = WFluxanlTool() #alpha == 0.01
        high_alpha_tester = WFluxanlTool() #alpha == 0.20

        out_low = self.mktemp(prefix='fluxanl')
        out_high = self.mktemp(prefix='fluxanl')

        default_alpha_tester.make_parser_and_process(args=['--alpha=0.05', '-o={}'.format(self.outfile)])
        low_alpha_tester.make_parser_and_process(args=['--alpha=0.01', '--autocorrel-alpha=0.05', '-o={}'.format(out_low)])
        high_alpha_tester.make_parser_and_process(args=['--alpha=0.20', '--autocorrel-alpha=0.05', '-o={}'.format(out_high)])

        assert default_alpha_tester.alpha == 0.05
        assert low_alpha_tester.alpha == 0.01
        assert high_alpha_tester.alpha == 0.20

        default_alpha_tester.go()
        low_alpha_tester.go()
        high_alpha_tester.go()

        f_default = h5py.File(self.outfile)
        f_low = h5py.File(out_low)
        f_high = h5py.File(out_high)

        data_default = get_flux_group(f_default)
        data_low = get_flux_group(f_low)
        data_high = get_flux_group(f_high)

        ci_len_default = data_default['mean_flux_ci_ub'][0] - data_default['mean_flux_ci_lb'][0]
        ci_len_low = data_low['mean_flux_ci_ub'][0] - data_low['mean_flux_ci_lb'][0]
        ci_len_high = data_high['mean_flux_ci_ub'][0] - data_high['mean_flux_ci_lb'][0]

        assert data_default['mean_flux'][0] == data_low['mean_flux'][0] == data_high['mean_flux'][0], 'Mean fluxes differ at different alphas'
        assert ci_len_low >= ci_len_default >= ci_len_high
