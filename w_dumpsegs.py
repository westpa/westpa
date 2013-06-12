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
from westtools import WESTTool, WESTDataReader
from west import Segment

import logging
log = logging.getLogger('westtools.w_dumpsegs')

class WDumpSegs(WESTTool):
    prog='w_dumpsegs'
    description = '''\
Dump segment data as text. This is very inefficient, so this tool should be used
as a last resort (use hdfview/h5ls to look at data, and access HDF5 directly for
significant analysis tasks). 
'''
    
    def __init__(self):
        super(WDumpSegs,self).__init__()
        self.data_reader = WESTDataReader()
        self.n_iter = None
        self.output_file = None
        self.print_pcoords = False
    
    def add_args(self, parser):
        self.data_reader.add_args(parser)


        parser.add_argument('-p', '--print-pcoords', dest='print_pcoords', action='store_true',
                            help='print initial and final progress coordinates for each segment')
        parser.add_argument('-i', '--iteration', dest='n_iter', type=int,
                            help='Use data from iteration N_ITER (default: last complete iteration)')
        parser.add_argument('-o', '--output', dest='output_file',
                            help='Store output in OUTPUT_FILE (default: write to standard output).')
        

    def process_args(self, args):
        self.data_reader.process_args(args)
        self.data_reader.open()
        self.n_iter = args.n_iter or self.data_reader.current_iteration-1
        self.output_file = open(args.output_file, 'wt') if args.output_file else sys.stdout
        self.print_pcoords = args.print_pcoords
    
    def go(self):
        segments = self.data_reader.get_segments(self.n_iter)
        
        max_seg_id_len = len(str(max(segment.seg_id for segment in segments)))
        max_status_name_len = max(map(len,Segment.status_names.itervalues()))
        max_endpoint_type_len = max(map(len,Segment.endpoint_type_names.itervalues()))
        max_n_parents_len = len(str(max(len(segment.wtg_parent_ids) for segment in segments)))
        
        report_line = ( '{segment.n_iter:d}  {segment.seg_id:{max_seg_id_len}d}  {segment.weight:20.14g}' 
                        +'  {status_name:{max_status_name_len}s} ({segment.status})'
                        +'  {segment.walltime:<12.6g} {segment.cputime:<12.6g}'
                        +'  {endpoint_type_name:{max_endpoint_type_len}s} ({segment.endpoint_type})'
                        +'  {n_parents:{max_n_parents_len}d} {segment.parent_id:{max_seg_id_len}d} {parents_str}' 
                        +'\n')
        pcoord_lines = ('  pcoord[0]  = {init_pcoord}\n  pcoord[-1] = {final_pcoord}'
                        +'\n')
        for (_seg_id, segment) in enumerate(segments):    
            parents_str = '['+', '.join(map(str,sorted(segment.wtg_parent_ids)))+']'
            init_pcoord_str = '[' + ', '.join('{pcval:<12.6g}'.format(pcval=float(pce)) for pce in segment.pcoord[0]) + ']'
            final_pcoord_str = '[' + ', '.join('{pcval:<12.6g}'.format(pcval=float(pce)) for pce in segment.pcoord[-1]) + ']'
            self.output_file.write(report_line.format(segment=segment, 
                                                 status_name = segment.status_names[segment.status],
                                                 endpoint_type_name = segment.endpoint_type_names[segment.endpoint_type],
                                                 parents_str = parents_str,
                                                 n_parents = len(segment.wtg_parent_ids),
                                                 max_seg_id_len=max_seg_id_len,
                                                 max_status_name_len=max_status_name_len,
                                                 max_endpoint_type_len=max_endpoint_type_len,
                                                 max_n_parents_len=max_n_parents_len))
            if self.print_pcoords:
                self.output_file.write(pcoord_lines.format(init_pcoord = init_pcoord_str,final_pcoord = final_pcoord_str))

if __name__ == '__main__':
    WDumpSegs().main()

