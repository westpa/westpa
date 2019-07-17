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



import logging, warnings
log = logging.getLogger(__name__)

import itertools, re

import numpy, h5py

import west, westpa
from oldtools.aframe import AnalysisMixin
from west import Segment
from oldtools.miscfn import parse_int_list

class WESTDataReaderMixin(AnalysisMixin):
    '''A mixin for analysis requiring access to the HDF5 files generated during a WEST run.'''

    def __init__(self):
        super(WESTDataReaderMixin,self).__init__()
        
        self.data_manager = None
        self.west_h5name = None
        
        # Whether pcoord caching is active
        self.__cache_pcoords = False        

        # Cached items
        self.__c_summary = None
        self.__c_iter_groups = dict()
        self.__c_seg_id_ranges = dict()
        self.__c_seg_indices = dict()
        self.__c_wtg_parent_arrays = dict()
        self.__c_parent_arrays = dict()
        self.__c_pcoord_arrays = dict()
        self.__c_pcoord_datasets = dict()
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upcall = super(WESTDataReaderMixin,self).add_args
            except AttributeError:
                pass
            else:
                upcall(parser)
        
        group = parser.add_argument_group('WEST input data options')
        group.add_argument('-W', '--west-data', dest='west_h5name', metavar='WEST_H5FILE',
                           help='''Take WEST data from WEST_H5FILE (default: read from the HDF5 file specified in west.cfg).''')

    def process_args(self, args, upcall = True):            
        if args.west_h5name:
            self.west_h5name = args.west_h5name
        else:
            westpa.rc.config.require(['west','data','west_data_file'])
            self.west_h5name = westpa.rc.config.get_path(['west','data','west_data_file']) 
        
        westpa.rc.pstatus("Using WEST data from '{}'".format(self.west_h5name))
        
        self.data_manager = westpa.rc.get_data_manager()
        self.data_manager.backing_file = self.west_h5name        
        self.data_manager.open_backing(mode='r')
        
        if upcall:
            try:
                upfunc = super(WESTDataReaderMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)

    def clear_run_cache(self):
        del self.__c_summary
        del self.__c_iter_groups, self.__c_seg_id_ranges, self.__c_seg_indices, self.__c_parent_arrays, self.__c_parent_arrays
        del self.__c_pcoord_arrays, self.__c_pcoord_datasets
        
        self.__c_summary = None
        self.__c_iter_groups = dict()
        self.__c_seg_id_ranges = dict()
        self.__c_seg_indices = dict()
        self.__c_parent_arrays = dict()
        self.__c_wtg_parent_arrays = dict()
        self.__c_pcoord_arrays = dict()
        self.__c_pcoord_datasets = dict()

    @property
    def cache_pcoords(self):
        '''Whether or not to cache progress coordinate data. While caching this data
        can significantly speed up some analysis operations, this requires
        copious RAM.
        
        Setting this to False when it was formerly True will release any cached data. 
        '''
        return self.__cache_pcoords

    @cache_pcoords.setter
    def cache_pcoords(self, cache):
        self.__cache_pcoords = cache
        
        if not cache:
            del self.__c_pcoord_arrays
            self.__c_pcoord_arrays = dict()
            
    def get_summary_table(self):
        if self.__c_summary is None:
            self.__c_summary = self.data_manager.we_h5file['/summary'][...]
        return self.__c_summary

    def get_iter_group(self, n_iter):
        '''Return the HDF5 group corresponding to ``n_iter``'''
        try:
            return self.__c_iter_groups[n_iter]
        except KeyError:
            iter_group = self.data_manager.get_iter_group(n_iter)
            return iter_group
    
    def get_segments(self, n_iter, include_pcoords = True):
        '''Return all segments present in iteration n_iter'''
        return self.get_segments_by_id(n_iter, self.get_seg_ids(n_iter, None), include_pcoords)
    
    def get_segments_by_id(self, n_iter, seg_ids, include_pcoords = True):
        '''Get segments from the data manager, employing caching where possible'''
        
        if len(seg_ids) == 0: return []
        
        seg_index  = self.get_seg_index(n_iter)
        all_wtg_parent_ids = self.get_wtg_parent_array(n_iter)
        
        segments = []

        if include_pcoords:
            pcoords = self.get_pcoords(n_iter, seg_ids)
        
        for (isegid, seg_id) in enumerate(seg_ids):
            row = seg_index[seg_id]
            parents_offset = row['wtg_offset']
            n_parents = row['wtg_n_parents']
            segment = Segment(seg_id = seg_id,
                              n_iter = n_iter,
                              status = row['status'],
                              endpoint_type = row['endpoint_type'],
                              walltime = row['walltime'],
                              cputime = row['cputime'],
                              weight = row['weight'],
                              )
            if include_pcoords:
                segment.pcoord = pcoords[isegid]

            parent_ids = all_wtg_parent_ids[parents_offset:parents_offset+n_parents]
            segment.wtg_parent_ids = {int(parent_id) for parent_id in parent_ids}
            segment.parent_id = int(parent_ids[0])
            segments.append(segment)

        return segments        
    
    def get_children(self, segment, include_pcoords=True):
        parents = self.get_parent_array(segment.n_iter+1)
        seg_ids = self.get_seg_ids(segment.n_iter+1, parents == segment.seg_id)
        return self.get_segments_by_id(segment.n_iter+1, seg_ids, include_pcoords)

    def get_seg_index(self, n_iter):
        try:
            return self.__c_seg_indices[n_iter]
        except KeyError:
            seg_index = self.__c_seg_indices[n_iter] = self.get_iter_group(n_iter)['seg_index'][...]
            return seg_index
        
    def get_wtg_parent_array(self, n_iter):
        try:
            return self.__c_wtg_parent_arrays[n_iter]
        except KeyError:
            parent_array = self.__c_wtg_parent_arrays[n_iter] = self.get_iter_group(n_iter)['wtgraph'][...]
            return parent_array
                
    def get_parent_array(self, n_iter):
        try:
            return self.__c_parent_arrays[n_iter]
        except KeyError:
            parent_array = self.get_seg_index(n_iter)['parent_id']
            self.__c_parent_arrays[n_iter] = parent_array
            return parent_array
                
    def get_pcoord_array(self, n_iter):
        try:
            return self.__c_pcoord_arrays[n_iter]
        except KeyError:
            pcoords = self.__c_pcoord_arrays[n_iter] = self.get_iter_group(n_iter)['pcoord'][...]
            return pcoords
        
    def get_pcoord_dataset(self, n_iter):
        try:
            return self.__c_pcoord_datasets[n_iter]
        except KeyError:
            pcoord_ds = self.__c_pcoord_datasets[n_iter] = self.get_iter_group(n_iter)['pcoord']
            return pcoord_ds
        
    def get_pcoords(self, n_iter, seg_ids):
        if self.__cache_pcoords:
            pcarray = self.get_pcoord_array(n_iter)
            return [pcarray[seg_id,...] for seg_id in seg_ids]
        else:
            return self.get_pcoord_dataset(n_iter)[list(seg_ids),...]
        
    def get_seg_ids(self, n_iter, bool_array = None):
        try:
            all_ids = self.__c_seg_id_ranges[n_iter]
        except KeyError:
            all_ids = self.__c_seg_id_ranges[n_iter] = numpy.arange(0,len(self.get_seg_index(n_iter)), dtype=numpy.uint32)
            
            
        if bool_array is None:             
            return all_ids
        else:        
            seg_ids = all_ids[bool_array]        
            try:
                if len(seg_ids) == 0: return []
            except TypeError:
                # Not iterable, for some bizarre reason
                return [seg_ids]
            else:
                return seg_ids
    
    def get_created_seg_ids(self, n_iter):
        '''Return a list of seg_ids corresponding to segments which were created for the given iteration (are not
        continuations).'''
        
        # Created segments have parent_id < 0
        parent_ids = self.get_parent_array(n_iter)        
        return self.get_seg_ids(n_iter, parent_ids < 0)

    def max_iter_segs_in_range(self, first_iter, last_iter):
        '''Return the maximum number of segments present in any iteration in the range selected'''
        n_particles = self.get_summary_table()['n_particles']
        return n_particles[first_iter-1:last_iter].max()
                
    def total_segs_in_range(self, first_iter, last_iter):
        '''Return the total number of segments present in all iterations in the range selected'''
        n_particles = self.get_summary_table()['n_particles']
        return n_particles[first_iter-1:last_iter].sum()

    def get_pcoord_len(self, n_iter):
        '''Get the length of the progress coordinate array for the given iteration.'''
        pcoord_ds = self.get_pcoord_dataset(n_iter)
        return pcoord_ds.shape[1]
    
    def get_total_time(self, first_iter = None, last_iter = None, dt=None):
        '''Return the total amount of simulation time spanned between first_iter and last_iter (inclusive).'''
        first_iter = first_iter or self.first_iter
        last_iter = last_iter or self.last_iter
        dt = dt or getattr(self, 'dt', 1.0)
        
        total_len = 0
        for n_iter in range(first_iter, last_iter+1):
            total_len += self.get_pcoord_len(n_iter) - 1
        return total_len * dt 
    

class ExtDataReaderMixin(AnalysisMixin):
    '''An external data reader, primarily designed for reading brute force data, but also suitable
    for any auxiliary datasets required for analysis.
    '''
    
    default_chunksize = 8192
    
    def __init__(self):
        super(ExtDataReaderMixin,self).__init__()
        
        self.ext_input_nargs = '+'
        self.ext_input_filenames = []
        self.ext_input_chunksize = self.default_chunksize
        self.ext_input_usecols = None
        self.ext_input_comment_regexp = None
        self.ext_input_sep_regexp = None
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upcall = super(ExtDataReaderMixin,self).add_args
            except AttributeError:
                pass
            else:
                upcall(parser)

        input_options = parser.add_argument_group('external data input options')        
        input_options.add_argument('datafiles', nargs=self.ext_input_nargs, metavar='DATAFILE',
                                   help='''Data file(s) to analyze, either text or Numpy (.npy or .npz) format.
                                   Uncompressed numpy files will be memory-mapped, allowing analysis of data larger than 
                                   available RAM (though not larger than the available address space).''')
        input_options.add_argument('--usecols', dest='usecols', metavar='COLUMNS', type=parse_int_list,
                                   help='''Use only the given COLUMNS from the input file(s), e.g. "0", "0,1", 
                                   "0:5,7,9:10".''')
        input_options.add_argument('--chunksize', dest='chunksize', type=int, default=self.default_chunksize,
                                   help='''Process input data in blocks of size CHUNKSIZE.  This will only reduce memory
                                   requirements when using uncompressed Numpy (.npy) format input. (Default: %(default)d.)''')    

    def process_args(self, args, upcall = True):            

        if args.usecols:
            westpa.rc.pstatus('Using only the following columns from external input: {!s}'.format(args.usecols))
            self.ext_input_usecols = args.usecols
        else:
            self.ext_input_usecols = None
            
        self.ext_input_filenames = args.datafiles
        self.ext_input_chunksize = args.chunksize or self.default_chunksize

        if upcall:
            try:
                upfunc = super(ExtDataReaderMixin,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)

    def is_npy(self, filename):
        with open(filename, 'rb') as fileobj:
            first_bytes = fileobj.read(len(numpy.lib.format.MAGIC_PREFIX))
        
        if first_bytes == numpy.lib.format.MAGIC_PREFIX:
            return True
        else:
            return False
                
    def load_npy_or_text(self, filename):
        '''Load an array from an existing .npy file, or read a text file and
        convert to a NumPy array.  In either case, return a NumPy array.  If a 
        pickled NumPy dataset is found, memory-map it read-only.  If the specified
        file does not contain a pickled NumPy array, attempt to read the file using
        numpy.loadtxt(filename).'''
        
        if self.is_npy(filename):
            return numpy.load(filename, 'r')
        else:
            return numpy.loadtxt(filename)
            
    def text_to_h5dataset(self, fileobj, group, dsname, dtype=numpy.float64, 
                             skiprows=0, usecols=None, 
                             chunksize=None):
        '''Read text-format data from the given filename or file-like object ``fileobj`` and write to a newly-created dataset
        called ``dsname`` in the HDF5 group ``group``.  The data is stored as type ``dtype``.  By default, the shape is
        taken as (number of lines, number of columns); columns can be omitted by specifying a list for ``usecols``,
        and lines can be skipped by using ``skiprows``.  Data is read in chunks of ``chunksize`` rows.'''
        
        try:
            fileobj.readline
        except AttributeError:
            fileobj = open(fileobj, 'rt')
            
        usecols = usecols or self.usecols
        chunksize = chunksize or self.ext_input_chunksize
        
        linenumber = 0
        for iskip in range(skiprows or 0):
            fileobj.readline()
            linenumber += 1
        
        nrows = 0
        irow = 0
        ncols_input = None # number of columns in input
        ncols_store = None # number of columns to store
        databuffer = None
        dataset = None
        
        re_split_comments = self.ext_input_comment_regexp
        re_split_fields = self.ext_input_sep_regexp
        
        for line in fileobj:
            linenumber += 1
            
            # Discard comments and extraneous whitespace
            if re_split_comments is not None:
                record_text = re_split_comments.split(line, 1)[0].strip()
            else:
                record_text = line.split('#', 1)[0].strip()
                
            if not record_text:
                continue
            
            if re_split_fields is not None:
                fields = re_split_fields.split(record_text)
            else:
                fields = record_text.split()
            
            # Check that the input size hasn't change (blank lines excluded)
            if not ncols_input:
                ncols_input = len(fields)
            elif len(fields) != ncols_input:
                raise ValueError('expected {:d} columns at line {:d}, but found {:d}'
                                 .format(ncols_input, linenumber, len(fields)))

            # If this is the first time through the loop, allocate temporary storage
            if not ncols_store:
                ncols_store = len(usecols)
                databuffer = numpy.empty((chunksize, ncols_store), dtype)
                dataset = group.create_dataset(dsname, 
                                               shape=(0,ncols_store), maxshape=(None,ncols_store), chunks=(chunksize,ncols_store), 
                                               dtype=dtype)
            
            if usecols:
                for (ifield,iifield) in enumerate(usecols):
                    databuffer[irow,ifield] = dtype(fields[iifield])
            else:
                for (ifield, field) in enumerate(fields):
                    databuffer[irow,ifield] = dtype(field)

            nrows+=1
            irow+=1
            
            # Flush to HDF5 if necessary
            if irow == chunksize:
                westpa.rc.pstatus('\r  Read {:d} rows'.format(nrows), end='')
                westpa.rc.pflush()
                dataset.resize((nrows, ncols_store))
                dataset[-irow:] = databuffer
                irow = 0

        # Flush last bit                
        if irow > 0:
            dataset.resize((nrows, ncols_store))
            dataset[-irow:] = databuffer[:irow]
        westpa.rc.pstatus('\r  Read {:d} rows'.format(nrows))
        westpa.rc.pflush()
            
    def npy_to_h5dataset(self, array, group, dsname, usecols=None, chunksize=None):
        '''Store the given array into a newly-created dataset named ``dsname`` in the HDF5 group
        ``group``, optionally only storing a subset of columns.  Data is written ``chunksize`` rows at a time,
        allowing very large memory-mapped arrays to be copied.'''
        
        usecols = usecols or self.ext_input_usecols
        chunksize = chunksize or self.ext_input_chunksize
        
        if usecols:
            shape = (len(array),) + array[0][usecols].shape[1:]
        else:
            shape = array.shape
        
        if len(shape) == 1:
            shape = shape + (1,)
        maxlen = len(array)
        mw = len(str(maxlen))
        dataset = group.create_dataset(dsname, shape=shape, dtype=array.dtype)
        
        if usecols:
            for istart in range(0,maxlen,chunksize):
                iend = min(istart+chunksize,maxlen)
                dataset[istart:iend] = array[istart:iend, usecols]
                westpa.rc.pstatus('\r  Read {:{mw}d}/{:>{mw}d} rows'.format(iend,maxlen, mw=mw), end='')
                westpa.rc.pflush()
        else:
            for istart in range(0,maxlen,chunksize):
                dataset[istart:iend] = array[istart:iend]
                westpa.rc.pstatus('\r  Read {:{mw}d}/{:>{mw}d} rows'.format(iend,maxlen, mw=mw), end='')
                westpa.rc.pflush()                
        westpa.rc.pstatus()
        
        
                
class BFDataManager(AnalysisMixin):
    '''A class to manage brute force trajectory data.  The primary purpose is to read in and 
    manage brute force progress coordinate data for one or more trajectories.  The trajectories need not
    be the same length, but they do need to have the same time spacing for progress coordinate values.'''
        
    traj_index_dtype = numpy.dtype( [ ('pcoord_len', numpy.uint64),
                                      ('source_data', h5py.special_dtype(vlen=str)) ] )
    
    def __init__(self):
        super(BFDataManager,self).__init__()
        self.bf_h5name = None
        self.bf_h5file = None
        
    def add_args(self, parser, upcall = True):
        if upcall:
            try:
                upcall = super(BFDataManager,self).add_args
            except AttributeError:
                pass
            else:
                upcall(parser)
        
        group = parser.add_argument_group('brute force input data options')
        group.add_argument('-B', '--bfdata', '--brute-force-data', dest='bf_h5name', metavar='BF_H5FILE', default='bf_system.h5',
                           help='''Brute force data is/will be stored in BF_H5FILE (default: %(default)s).''')

    def process_args(self, args, upcall = True):
        self.bf_h5name = args.bf_h5name
        westpa.rc.pstatus("Using brute force data from '{}'".format(self.bf_h5name))
        
        if upcall:
            try:
                upfunc = super(BFDataManager,self).process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
        
    def _get_traj_group_name(self, traj_id):
        return 'traj_{:09d}'.format(traj_id)
    
    def update_traj_index(self, traj_id, pcoord_len, source_data):
        self.bf_h5file['traj_index'][traj_id] = (pcoord_len, source_data)
                
    def get_traj_group(self, traj_id):
        return self.bf_h5file[self._get_traj_group_name(traj_id)]
    
    def create_traj_group(self):
        new_traj_id = self.get_n_trajs()
        group = self.bf_h5file.create_group(self._get_traj_group_name(new_traj_id))
        self.bf_h5file['traj_index'].resize((new_traj_id+1,))
        return (new_traj_id, group)
        
    def get_n_trajs(self):
        return self.bf_h5file['traj_index'].shape[0]
    
    def get_traj_len(self,traj_id):
        return self.bf_h5file['traj_index'][traj_id]['pcoord_len']
    
    def get_max_traj_len(self):
        return self.bf_h5file['traj_index']['pcoord_len'].max()

    def get_pcoord_array(self, traj_id):
        return self.get_traj_group(traj_id)['pcoord'][...]
        
    def get_pcoord_dataset(self, traj_id):
        return self.get_traj_group(traj_id)['pcoord']
    
    def require_bf_h5file(self):
        if self.bf_h5file is None:
            assert self.bf_h5name
            self.bf_h5file = h5py.File(self.bf_h5name)
            try:
                self.bf_h5file['traj_index']
            except KeyError:
                # A new file; create the trajectory index
                self.bf_h5file.create_dataset('traj_index', shape=(0,), maxshape=(None,), dtype=self.traj_index_dtype)
        return self.bf_h5file
    
    def close_bf_h5file(self):
        if self.bf_h5file is not None:
            self.bf_h5file.close()
            self.bf_h5file = None
            

        
