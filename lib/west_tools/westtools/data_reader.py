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

import numpy
from numpy import index_exp
from .core import WESTToolComponent
import westpa
from westpa.extloader import get_object
from westpa.h5io import FnDSSpec, MultiDSSpec, SingleSegmentDSSpec, SingleIterDSSpec


def _get_parent_ids(n_iter, iter_group):
    seg_index = iter_group['seg_index']
    try:
        return seg_index['parent_id'][:]
    except ValueError:
        # field not found
        offsets = seg_index['parents_offset'][:]
        all_parents = iter_group['parents'][...]
        return numpy.require(all_parents.take(offsets),dtype=numpy.int64)
    else:
        return seg_index['parent_id']
    

class WESTDataReader(WESTToolComponent):
    '''Tool for reading data from WEST-related HDF5 files. Coordinates finding
    the main HDF5 file from west.cfg or command line arguments, caching of certain
    kinds of data (eventually), and retrieving auxiliary data sets from various
    places.'''
    
    def __init__(self):
        super(WESTDataReader,self).__init__()
        self.data_manager = westpa.rc.get_data_manager() 
        self.we_h5filename = None
        
        self._weight_dsspec = None
        self._parent_id_dsspec = None
        
    def add_args(self, parser):
        group = parser.add_argument_group('WEST input data options')
        group.add_argument('-W', '--west-data', dest='we_h5filename', metavar='WEST_H5FILE',
                           help='''Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in west.cfg).''')
        
    def process_args(self, args):
        if args.we_h5filename:
            self.data_manager.we_h5filename = self.we_h5filename = args.we_h5filename
        else:
            self.we_h5filename = self.data_manager.we_h5filename
        
    def open(self, mode='r'):
        self.data_manager.open_backing(mode)
        
    def close(self):
        self.data_manager.close_backing()
        
    def __getattr__(self, key):
        return getattr(self.data_manager, key)
    
    @property
    def weight_dsspec(self):
        if self._weight_dsspec is None:
            assert self.we_h5filename is not None
            self._weight_dsspec = SingleIterDSSpec(self.we_h5filename, 'seg_index', slice=index_exp['weight'])
        return self._weight_dsspec

    @property
    def parent_id_dsspec(self):
        if self._parent_id_dsspec is None:
            assert self.we_h5filename is not None
            #self._parent_id_dsspec = SingleIterDSSpec(self.we_h5filename, 'seg_index', slice=index_exp['parent_id'])
            self._parent_id_dsspec = FnDSSpec(self.we_h5filename, _get_parent_ids)
        return self._parent_id_dsspec

    def __enter__(self):
        self.open('r')
        return self
        
    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.close()
        return False
            

class WESTDSSynthesizer(WESTToolComponent):
    '''Tool for synthesizing a dataset for analysis from other datasets. This
    may be done using a custom function, or a list of "data set specifications".
    It is anticipated that if several source datasets are required, then a tool
    will have multiple instances of this class.'''
    
    group_name = 'input dataset options'
    
    def __init__(self, default_dsname = None, h5filename=None):
        super(WESTDSSynthesizer,self).__init__()

        self.h5filename = h5filename
        self.default_dsname = default_dsname
        
        self.dsspec = None
        
    def add_args(self, parser):
        igroup = parser.add_argument_group(self.group_name).add_mutually_exclusive_group(required=not bool(self.default_dsname))

        igroup.add_argument('--construct-dataset',
                            help='''Use the given function (as in module.function) to extract source data.
                            This function will be called once per iteration as function(n_iter, iter_group)
                            to construct data for one iteration. Data returned must be indexable as
                            [seg_id][timepoint][dimension]''')
        
        igroup.add_argument('--dsspecs', nargs='+', metavar='DSSPEC',
                            help='''Construct source data from one or more DSSPECs.''')
        
    def process_args(self, args):
        if args.construct_dataset:
            self.dsspec = FnDSSpec(self.h5filename, get_object(args.construct_dataset,path=['.']))
        elif args.dsspecs:
            self.dsspec = MultiDSSpec([SingleSegmentDSSpec.from_string(dsspec, self.h5filename)
                                       for dsspec in args.dsspecs])
        else:
            # we can only get here if a default dataset name was specified
            assert self.default_dsname
            self.dsspec = SingleSegmentDSSpec(self.h5filename, self.default_dsname)
        