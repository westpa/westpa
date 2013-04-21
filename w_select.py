from __future__ import print_function, division; __metaclass__ = type
import sys
from westtools.tool_classes import WESTTool, WESTDataReader
from itertools import imap
import numpy

import westpa
from west.data_manager import (weight_dtype, n_iter_dtype, seg_id_dtype, utime_dtype, vstr_dtype, 
                               istate_type_dtype, istate_status_dtype)
from westpa.extloader import get_object

def all_segments(n_iter, iter_group):
    return numpy.arange(0, iter_group['seg_index'].shape[0], 1, dtype=seg_id_dtype)
    
def find_matching_segments(predicate, invert_predicate=False, ancestors=True, data_manager=None):
    data_manager = data_manager or westpa.rc.get_data_manager()
    
    match_pairs = []    
    parent_ids = set()

    for n_iter in xrange(data_manager.current_iteration-1, 0, -1):
        iter_group = data_manager.get_iter_group(n_iter)
        
        matching_ids = set(imap(long, predicate(n_iter, iter_group)))

        if invert_predicate:
            all_ids = set(imap(long,xrange(0, iter_group['seg_index'].shape[0])))
            matching_ids = all_ids - matching_ids
        
        matching_ids.update(parent_ids)
                    
        if matching_ids:                
            if ancestors:
                parent_ids = set(imap(long,data_manager.get_parent_ids(n_iter, matching_ids)))
            else:
                parent_ids = set()
            
            match_pairs.extend((n_iter, seg_id) for seg_id in matching_ids)
            
        del iter_group
            
    return match_pairs
        
    
class WSelectTool(WESTTool):
    prog='w_select'
    description = '''\
Select dynamics segments matching various criteria, such as "all segments
in trajectories ending in recycling". The predicate function used to select 
segments may be from a pre-defined selection, or user-specified via a 
Python callable.  Ancestors of matching segments may be included or
omitted.  The list of matching segments is emitted as a sequence of
N_ITER:SEG_ID pairs, separated by newlines.
'''
    
    def __init__(self):
        super(WSelectTool,self).__init__()
        
        self.data_reader = WESTDataReader()
        self.output_file = None
        self.predicate = None
        self.invert = False
        self.ancestors = True

        
    # Interface for command-line tools
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        
        sgroup = parser.add_argument_group('selection options')
        fngroup = sgroup.add_mutually_exclusive_group()
        fngroup.add_argument('-p', '--predicate-function', metavar='MODULE.FUNCTION',
                             help='''Use the given predicate function to match segments. This function
                             should take an iteration number and the HDF5 group corresponding to that
                             iteration and return a sequence of seg_ids matching the predicate, as in
                             ``match_predicate(n_iter, iter_group)``.''')
        fngroup.add_argument('--recycled', nargs='?', metavar='TSTATE_INDEX',
                             help='''Match segments which are recycled. If the optional TSTATE_INDEX is
                             given, match only segments recycled from the given target state. This is the
                             predicate to use when searching for segments involved in "successful"
                             trajectories.''')
        
        sgroup.add_argument('-v', '--invert', dest='invert', action='store_true',
                            help='''Invert the match predicate.''')
        sgroup.add_argument('--no-ancestors', action='store_true',
                            help='''Consider only segments matching the predicate function, not their 
                            ancestors. (By default, segments matching the predicate and their ancestors
                            are output.)''')

        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('-o', '--output',
                            help='''Write output to OUTPUT (default: standard output).''')
    
    def process_args(self, args): 
        if args.predicate_function:
            predicate = get_object(args.predicate_function,path=['.'])
            if not callable(predicate):
                raise TypeError('predicate object {!r} is not callable'.format(predicate))
        else:
            predicate = all_segments
        self.predicate = predicate  
               
        self.invert = bool(args.invert)
        
        if args.no_ancestors:
            self.ancestors = False
        else:
            self.ancestors = True   

        if args.output is not None:
            self.output_file = open(args.output,'wt')
        else:
            self.output_file = sys.stdout
        
        self.data_reader.process_args(args)
        
    def go(self):
        self.data_reader.open('r')
        
        matching_pairs = find_matching_segments(self.predicate, invert_predicate=self.invert, 
                                                ancestors=self.ancestors, data_manager=self.data_reader.data_manager)
    
        output_file = self.output_file
        for (n_iter,seg_id) in matching_pairs:
            output_file.write('{:d}:{:d}\n'.format(n_iter,seg_id))
            
if __name__ == '__main__':
    WSelectTool().main()
    
