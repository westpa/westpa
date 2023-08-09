import itertools
import logging
import re

import numpy as np

import westpa
from westpa.oldtools.aframe import AnalysisMixin

log = logging.getLogger(__name__)


class KineticsAnalysisMixin(AnalysisMixin):
    def __init__(self):
        super().__init__()

        self.dt = None
        self.analysis_initial_bins = None
        self.analysis_final_bins = None

    def add_args(self, parser, upcall=True):
        if upcall:
            try:
                upfunc = super().add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)

        group = parser.add_argument_group('kinetics analysis options')
        group.add_argument(
            '--dt', dest='dt', type=float, default=1.0, help='Assume input data has a time spacing of DT (default: %(default)s).'
        )
        group.add_argument(
            '--initial-bins',
            dest='ibins_string',
            metavar='ILIST',
            help='''Only calculate statistics for transitions starting in bin ILIST.  This may be specified as a
                             comma-separated list of integers or ranges, as in "0,2-4,5,9"''',
        )
        group.add_argument(
            '--final-bins',
            dest='fbins_string',
            metavar='FLIST',
            help='''Only calculate statistics for transitions ending in bin FLIST.  This may be specified as a
                             comma-separated list of integers or ranges, as in "0,2-4,5,9"''',
        )

    def process_args(self, args, upcall=True):
        self.dt = args.dt
        westpa.rc.pstatus('Assuming input data timestep of {:g}'.format(self.dt))

        if args.ibins_string:
            self.analysis_initial_bins = self.parse_bin_range(args.ibins_string)
            westpa.rc.pstatus(
                'Will calculate kinetics data from transitions beginning in the following bins: {!s}'.format(
                    sorted(self.analysis_initial_bins)
                )
            )
        else:
            westpa.rc.pstatus('Will calculate kinetics data from transitions beginning in any bin.')

        if args.fbins_string:
            self.analysis_final_bins = self.parse_bin_range(args.fbins_string)
            westpa.rc.pstatus(
                'Will calculate kinetics data from transitions ending in the following bins: {!s}'.format(
                    sorted(self.analysis_final_bins)
                )
            )
        else:
            westpa.rc.pstatus('Will calculate kinetics data from transitions ending in any bin.')

        if upcall:
            try:
                upfunc = super().process_args
            except AttributeError:
                pass
            else:
                upfunc(args)

    def parse_bin_range(self, range_string):
        try:
            entries = set()
            fields = re.split(r'\s*,\s*', range_string)
            for field in fields:
                if '-' in field:
                    lb, ub = list(map(int, re.split(r'\s*-\s*', field)))
                    entries.update(list(range(lb, ub + 1)))
                else:
                    entries.add(int(field))
        except (ValueError, TypeError):
            raise ValueError('invalid bin range string {!r}'.format(range_string))
        else:
            return entries

    def check_bin_selection(self, n_bins=None):
        '''Check to see that the bin ranges selected by the user conform to the available bins (i.e.,
        bin indices are within the permissible range).  Also assigns the complete bin range if the
        user has not explicitly limited the bins to be considered.'''

        n_bins = n_bins or self.n_bins

        if self.analysis_initial_bins:
            if (np.array(list(self.analysis_initial_bins)) >= n_bins).any():
                raise ValueError('One or more initial bin indices is out of range.')
        else:
            self.analysis_initial_bins = set(range(n_bins))

        if self.analysis_final_bins:
            if (np.array(list(self.analysis_final_bins)) >= n_bins).any():
                raise ValueError('One or more final bin indices is out of range.')
        else:
            self.analysis_final_bins = set(range(n_bins))

    @property
    def selected_bin_pair_iter(self):
        return (tuple(pair) for pair in itertools.product(self.analysis_initial_bins, self.analysis_final_bins))
