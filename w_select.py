# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division; __metaclass__ = type
import sys
from westtools import WESTParallelTool, WESTDataReader, IterRangeSelection, ProgressIndicatorComponent
from itertools import imap
import numpy

import westpa
from westpa import h5io
from west.data_manager import seg_id_dtype, n_iter_dtype, weight_dtype
from westpa.extloader import get_object

def _find_matching_segments(west_datafile_name, n_iter, predicate, invert=False):
    '''Find all segments in iteration ``n_iter`` that match (or do not match, if
    ``invert`` is true) the given ``predicate``. Returns a sequence of matching
    seg_ids.'''

    with h5io.WESTPAH5File(west_datafile_name, 'r') as west_datafile:
        iter_group = west_datafile.get_iter_group(n_iter)
        nsegs = iter_group['seg_index'].shape[0]
        matching_ids = set(imap(long, predicate(n_iter, iter_group)))

        if invert:
            matching_ids = set(xrange(nsegs)) - matching_ids

        matchvec = numpy.fromiter(matching_ids, dtype=seg_id_dtype, count=len(matching_ids))
        matchvec.sort()
        return n_iter, matchvec


class WSelectTool(WESTParallelTool):
    prog='w_select'
    description = '''\
Select dynamics segments matching various criteria. This requires a
user-provided prediate function. By default, only matching segments are
stored. If the -a/--include-ancestors option is given, then matching segments
and their ancestors will be stored.


-----------------------------------------------------------------------------
Predicate function
-----------------------------------------------------------------------------

Segments are selected based on a predicate function, which must be callable
as ``predicate(n_iter, iter_group)`` and return a collection of segment IDs
matching the predicate in that iteration.

The predicate may be inverted by specifying the -v/--invert command-line
argument.


-----------------------------------------------------------------------------
Output format
-----------------------------------------------------------------------------

The output file (-o/--output, by default "select.h5") contains the following
datasets:

  ``/n_iter`` [iteration]
    *(Integer)* Iteration numbers for each entry in other datasets.

  ``/n_segs`` [iteration]
    *(Integer)* Number of segment IDs matching the predicate (or inverted
    predicate, if -v/--invert is specified) in the given iteration.

  ``/seg_ids`` [iteration][segment]
    *(Integer)* Matching segments in each iteration. For an iteration
    ``n_iter``, only the first ``n_iter`` entries are valid. For example,
    the full list of matching seg_ids in the first stored iteration is
    ``seg_ids[0][:n_segs[0]]``.

  ``/weights`` [iteration][segment]
    *(Floating-point)* Weights for each matching segment in ``/seg_ids``.


-----------------------------------------------------------------------------
Command-line arguments
-----------------------------------------------------------------------------
'''

    def __init__(self):
        super(WSelectTool,self).__init__()

        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection()
        self.progress = ProgressIndicatorComponent()
        self.output_file = None
        self.output_filename = None
        self.predicate = None
        self.invert = False
        self.include_ancestors = False

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        sgroup = parser.add_argument_group('selection options')
        sgroup.add_argument('-p', '--predicate-function', metavar='MODULE.FUNCTION',
                             help='''Use the given predicate function to match segments. This function
                             should take an iteration number and the HDF5 group corresponding to that
                             iteration and return a sequence of seg_ids matching the predicate, as in
                             ``match_predicate(n_iter, iter_group)``.''')
        sgroup.add_argument('-v', '--invert', dest='invert', action='store_true',
                            help='''Invert the match predicate.''')
        sgroup.add_argument('-a', '--include-ancestors', action ='store_true',
                            help='''Include ancestors of matched segments in output.''')

        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output', default='select.h5',
                            help='''Write output to OUTPUT (default: %(default)s).''')
        self.progress.add_args(parser)

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)

        predicate = get_object(args.predicate_function,path=['.'])
        if not callable(predicate):
            raise TypeError('predicate object {!r} is not callable'.format(predicate))
        self.predicate = predicate
        self.invert = bool(args.invert)
        self.include_ancestors = bool(args.include_ancestors)
        self.output_filename = args.output

    def go(self):
        self.data_reader.open('r')
        output_file = h5io.WESTPAH5File(self.output_filename, mode='w')
        pi = self.progress.indicator

        iter_start, iter_stop = self.iter_range.iter_start, self.iter_range.iter_stop
        iter_count = iter_stop - iter_start

        output_file.create_dataset('n_iter', dtype=n_iter_dtype, data=range(iter_start,iter_stop))
        current_seg_count = 0
        seg_count_ds = output_file.create_dataset('n_segs', dtype=numpy.uint, shape=(iter_count,))
        matching_segs_ds = output_file.create_dataset('seg_ids', shape=(iter_count,0), maxshape=(iter_count,None),
                                                      dtype=seg_id_dtype,
                                                      chunks=h5io.calc_chunksize((iter_count,1000000), seg_id_dtype),
                                                      shuffle=True, compression=9)
        weights_ds = output_file.create_dataset('weights', shape=(iter_count,0), maxshape=(iter_count,None),
                                                dtype=weight_dtype,
                                                chunks=h5io.calc_chunksize((iter_count,1000000), weight_dtype),
                                                shuffle=True,compression=9)

        with pi:
            pi.new_operation('Finding matching segments', extent=iter_count)
#             futures = set()
#             for n_iter in xrange(iter_start,iter_stop):
#                 futures.add(self.work_manager.submit(_find_matching_segments, 
#                                                      args=(self.data_reader.we_h5filename,n_iter,self.predicate,self.invert)))

#             for future in self.work_manager.as_completed(futures):
            for future in self.work_manager.submit_as_completed(((_find_matching_segments,
                                                                  (self.data_reader.we_h5filename,n_iter,self.predicate,self.invert),
                                                                  {}) for n_iter in xrange(iter_start,iter_stop)),
                                                                self.max_queue_len):
                n_iter, matching_ids = future.get_result()
                n_matches = len(matching_ids)

                if n_matches:
                    if n_matches > current_seg_count:
                        current_seg_count = len(matching_ids)
                        matching_segs_ds.resize((iter_count,n_matches))
                        weights_ds.resize((iter_count,n_matches))
                        current_seg_count = n_matches

                    seg_count_ds[n_iter-iter_start] = n_matches
                    matching_segs_ds[n_iter-iter_start,:n_matches] = matching_ids
                    weights_ds[n_iter-iter_start,:n_matches] = self.data_reader.get_iter_group(n_iter)['seg_index']['weight'][sorted(matching_ids)]
                del matching_ids
                pi.progress += 1

            if self.include_ancestors:
                pi.new_operation('Tracing ancestors of matching segments', extent=iter_count)
                from_previous = set()
                current_seg_count = matching_segs_ds.shape[1]
                for n_iter in xrange(iter_stop-1, iter_start-1, -1):
                    iiter = n_iter - iter_start
                    n_matches = seg_count_ds[iiter]
                    matching_ids = set(from_previous)
                    if n_matches:
                        matching_ids.update(matching_segs_ds[iiter, :seg_count_ds[iiter]])
                    from_previous.clear()

                    n_matches = len(matching_ids)
                    if n_matches > current_seg_count:
                        matching_segs_ds.resize((iter_count,n_matches))
                        weights_ds.resize((iter_count,n_matches))
                        current_seg_count = n_matches

                    if n_matches > 0:
                        seg_count_ds[iiter] = n_matches
                        matching_ids = sorted(matching_ids)
                        matching_segs_ds[iiter,:n_matches] = matching_ids
                        weights_ds[iiter,:n_matches] = self.data_reader.get_iter_group(n_iter)['seg_index']['weight'][sorted(matching_ids)]
                        parent_ids = self.data_reader.get_iter_group(n_iter)['seg_index']['parent_id'][sorted(matching_ids)]
                        from_previous.update(parent_id for parent_id in parent_ids if parent_id >= 0) # filter initial states
                        del parent_ids
                    del matching_ids
                    pi.progress += 1

if __name__ == '__main__':
    WSelectTool().main()
    
