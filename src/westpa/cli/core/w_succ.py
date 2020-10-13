import argparse
import sys

import numpy as np

import westpa

from westpa.core.segment import Segment
from westpa.oldtools.aframe import WESTAnalysisTool, WESTDataReaderMixin, CommonOutputMixin

import logging

log = logging.getLogger('w_succ')


class WSucc(CommonOutputMixin, WESTDataReaderMixin, WESTAnalysisTool):
    def __init__(self):
        super().__init__()
        self.include_args['CommonOutputMixin']['print_bin_labels'] = False
        self.output_file = sys.stdout

    def find_successful_trajs(self):
        pcoord_formats = {
            'u8': '%20d',
            'i8': '%20d',
            'u4': '%10d',
            'i4': '%11d',
            'u2': '%5d',
            'i2': '%6d',
            'f4': '%14.7g',
            'f8': '%23.15g',
        }

        if not self.output_suppress_headers:
            self.output_file.write(
                '''\
# successful (recycled) segments
# column 0:    iteration
# column 1:    seg_id
# column 2:    weight
# column>2:    final progress coordinate value
'''
            )
        for n_iter in range(1, self.data_manager.current_iteration):
            seg_index = self.get_seg_index(n_iter)
            all_seg_ids = np.arange(len(seg_index), dtype=np.int_)
            recycled_seg_ids = all_seg_ids[seg_index[:]['endpoint_type'] == Segment.SEG_ENDPOINT_RECYCLED]

            if len(recycled_seg_ids) == 0:
                # Attemping to retrieve a 0-length selection from HDF5 (the pcoords below) fails
                continue

            pcoord_ds = self.get_pcoord_dataset(n_iter)
            pcoord_len = pcoord_ds.shape[1]
            pcoord_ndim = pcoord_ds.shape[2]
            final_pcoords = self.get_pcoord_dataset(n_iter)[recycled_seg_ids, pcoord_len - 1, :]
            # The above HDF5 selection always returns a vector; we want a 2-d array
            final_pcoords.shape = (len(recycled_seg_ids), pcoord_ndim)

            for (ipc, seg_id) in enumerate(recycled_seg_ids):
                self.output_file.write('%8d    %8d    %20.14g' % (n_iter, seg_id, seg_index[seg_id]['weight']))
                fields = ['']
                for field in final_pcoords[ipc]:
                    fields.append(pcoord_formats.get(field.dtype.str[1:], '%s') % field)
                self.output_file.write('    '.join(fields))
                self.output_file.write('\n')


def entry_point():
    wsucc = WSucc()

    parser = argparse.ArgumentParser(
        'w_succ',
        description='''\
    List segments which successfully reach a target state''',
    )

    westpa.rc.add_args(parser)
    wsucc.add_args(parser)

    parser.add_argument(
        '-o',
        '--output',
        dest='output_file',
        help='Store output in OUTPUT_FILE (default: write to standard output).',
        type=argparse.FileType('wt'),
        default=sys.stdout,
    )

    args = parser.parse_args()
    westpa.rc.process_args(args, config_required=False)
    wsucc.process_args(args)
    wsucc.output_file = args.output_file

    wsucc.find_successful_trajs()


if __name__ == '__main__':
    entry_point()
