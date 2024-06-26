import math
import pytest
import h5py
import numpy

from westpa.cli.tools.w_fluxanl import WFluxanlTool


def get_flux_group(f):
    '''f is an h5py file type object'''

    return f['target_flux']['index']


class Test_W_FLUXANL_NEW:

    mean_flux = None

    @pytest.mark.parametrize('evol_step', [None])
    @pytest.mark.parametrize('evol', [None, True])
    @pytest.mark.parametrize('nsets', [None, 500])
    @pytest.mark.parametrize('autocorrel_alpha', [None, 0.01])
    @pytest.mark.parametrize('alpha', [None, 0.01, 0.20])
    def test_args_alpha_autocorrel_nsets_evol_step(self, ref_50iter, alpha, autocorrel_alpha, nsets, evol, evol_step):
        '''Testing arg combos: w_fluxanl runs as expected'''
        outfile = "fluxanl.h5"

        args = {
            'alpha': alpha,
            'autocorrel-alpha': autocorrel_alpha,
            'nsets': nsets,
            'evol': evol,
            'evol-step': evol_step,
            'output': outfile,
        }

        arglist = [
            '--{}={}'.format(key, value) if not isinstance(value, bool) else '--{}'.format(key)
            for (key, value) in list(args.items())
            if value is not None
        ]

        w = WFluxanlTool()
        w.make_parser_and_process(args=arglist)

        # test_outcome, err = self.check_args_processed(**args_dict)
        w.go()
        self.check_output(outfile, **args)

    def check_output(self, outfile, **kwargs):
        '''Check that output file for w_fluxanl is set up correctly'''

        with h5py.File(outfile) as f:
            assert 'target_flux' in list(f.keys()), "'target flux' group not in output file"
            target_group = f['target_flux']
            assert 'index' in list(target_group.keys()), "'index' group not in output file"
            index_group = target_group['index']

            mean_flux = index_group['mean_flux'][0]

            if self.mean_flux is not None:
                assert mean_flux == self.mean_flux, "Actual mean flux ({}) is different than expected ({})".format(
                    mean_flux, self.mean_flux
                )
            else:
                self.mean_flux = mean_flux

            expected_alpha = kwargs['alpha'] or 0.05
            actual_alpha = index_group.attrs['mcbs_alpha']

            assert expected_alpha == actual_alpha, "Actual alpha ({}) does not match expected ({})".format(
                actual_alpha, expected_alpha
            )

            expected_acalpha = kwargs['autocorrel-alpha'] or expected_alpha
            actual_acalpha = index_group.attrs['mcbs_autocorrel_alpha']

            assert expected_acalpha == actual_acalpha, "Actual autocorrel alpha ({}) does not match expected ({})".format(
                actual_alpha, expected_alpha
            )

            expected_nsets = kwargs['nsets'] or 10 ** (math.ceil(-math.log10(expected_alpha)) + 1)
            actual_nsets = index_group.attrs['mcbs_n_sets']

            assert expected_nsets == actual_nsets, "Actual n_sets ({}) does not match expected ({})".format(
                actual_nsets, expected_nsets
            )


class Test_W_Fluxanl_System:
    '''System level tests for w_fluxanl'''

    def test_alpha(self, ref_50iter):
        '''Confidence interval size decreases as alpha increases'''

        self.w = WFluxanlTool()
        self.outfile = 'fluxanl.h5'

        default_alpha_tester = self.w  # alpha == 0.05
        low_alpha_tester = WFluxanlTool()  # alpha == 0.01
        high_alpha_tester = WFluxanlTool()  # alpha == 0.20

        out_low = 'fluxanl_low.h5'
        out_high = 'fluxanl_high.h5'

        default_alpha_tester.make_parser_and_process(args=['--alpha=0.05', '-o={}'.format(self.outfile)])
        low_alpha_tester.make_parser_and_process(args=['--alpha=0.01', '--autocorrel-alpha=0.05', '-o={}'.format(out_low)])
        high_alpha_tester.make_parser_and_process(args=['--alpha=0.20', '--autocorrel-alpha=0.05', '-o={}'.format(out_high)])

        assert default_alpha_tester.alpha == 0.05
        assert low_alpha_tester.alpha == 0.01
        assert high_alpha_tester.alpha == 0.20

        default_alpha_tester.go()
        low_alpha_tester.go()
        high_alpha_tester.go()

        with h5py.File(self.outfile) as f_default, h5py.File(out_low) as f_low, h5py.File(out_high) as f_high:
            data_default = get_flux_group(f_default)
            data_low = get_flux_group(f_low)
            data_high = get_flux_group(f_high)

            ci_len_default = data_default['mean_flux_ci_ub'][0] - data_default['mean_flux_ci_lb'][0]
            ci_len_low = data_low['mean_flux_ci_ub'][0] - data_low['mean_flux_ci_lb'][0]
            ci_len_high = data_high['mean_flux_ci_ub'][0] - data_high['mean_flux_ci_lb'][0]

            error_message = 'Mean fluxes differ at different alphas'

            assert numpy.isclose(data_default['mean_flux'][0], data_low['mean_flux'][0]), error_message

            assert numpy.isclose(data_default['mean_flux'][0], data_high['mean_flux'][0]), error_message

            assert numpy.isclose(data_low['mean_flux'][0], data_high['mean_flux'][0]), error_message

            assert ci_len_low >= ci_len_default >= ci_len_high
