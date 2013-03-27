'''Miscellaneous routines to help with HDF5 input and output of WEST-related data.'''

import sys, os, getpass, socket, time
import numpy, h5py
from numpy import index_exp

#
# Constants and globals
#
default_iter_prec=8

#
# Helper functions
# 
def calc_chunksize(shape, dtype, max_chunksize=262144):
    '''Calculate a chunk size for HDF5 data, anticipating that access will slice
    along lower dimensions sooner than higher dimensions.'''
        
    chunk_shape = list(shape)
    dtype = numpy.dtype(dtype)
    for idim in xrange(len(shape)):
        chunk_nbytes = numpy.multiply.reduce(chunk_shape)*dtype.itemsize
        while chunk_shape[idim] > 1 and chunk_nbytes > max_chunksize:
            chunk_shape[idim] >>= 1 # divide by 2
            chunk_nbytes = numpy.multiply.reduce(chunk_shape)*dtype.itemsize
            
        if chunk_nbytes <= max_chunksize:
            break

    chunk_shape = tuple(chunk_shape)
    return chunk_shape

#
# Group and datset manipulation functions
#

def create_hdf5_group(parent_group, groupname, replace=False, creating_program=None):
    '''Create (or delete and recreate) and HDF5 group named ``groupname`` within
    the enclosing Group (object) ``parent_group``. If ``replace`` is True, then
    the group is replaced if present; if False, then an error is raised if the
    group is present. After the group is created, HDF5 attributes are set using
    `stamp_creator_data`.
    '''
    
    if replace:
        try:
            del parent_group[groupname]
        except KeyError:
            pass
        
    newgroup = parent_group.create_group(groupname)
    stamp_creator_data(newgroup)
    return newgroup

#
# Group and dataset labeling functions
#

def stamp_creator_data(h5group, creating_program = None):
    '''Mark the following on the HDF5 group ``h5group``:
    
      :creation_program:   The name of the program that created the group
      :creation_user:      The username of the user who created the group
      :creation_hostname:  The hostname of the machine on which the group was created
      :creation_time:      The date and time at which the group was created, in the
                           current locale.
      :creation_unix_time: The Unix time (seconds from the epoch, UTC) at which the
                           group was created.
      
    This is meant to facilitate tracking the flow of data, but should not be considered
    a secure paper trail (after all, anyone with write access to the HDF5 file can modify
    these attributes).    
    '''
    now = time.time()
    attrs = h5group.attrs
    
    attrs['creation_program'] = creating_program or sys.argv[0] or 'unknown program'
    attrs['creation_user'] = getpass.getuser()
    attrs['creation_hostname'] = socket.gethostname()
    attrs['creation_unix_time'] = now
    attrs['creation_time'] = time.strftime('%c', time.localtime(now))

def get_creator_data(h5group):
    '''Read back creator data as written by ``stamp_creator_data``, returning a dictionary with
    keys as described for ``stamp_creator_data``. Missing fields are denoted with None. 
    The ``creation_time`` field is returned as a string.'''
    attrs = h5group.attrs
    d = dict()
    for attr in ['creation_program', 'creation_user', 'creation_hostname', 'creation_unix_time', 'creation_time']:
        d[attr] = attrs.get(attr)
    return d


###
# Iteration range metadata
###
def stamp_iter_range(h5object, start_iter, stop_iter):
    '''Mark that the HDF5 object ``h5object`` (dataset or group) contains data from iterations
    start_iter <= n_iter < stop_iter.'''
    h5object.attrs['iter_start'] = start_iter
    h5object.attrs['iter_stop']  = stop_iter
    
def get_iter_range(h5object):
    '''Read back iteration range data written by ``stamp_iter_range``'''
    return int(h5object.attrs['iter_start']), int(h5object.attrs['iter_stop'])

def stamp_iter_step(h5group, iter_step):
    '''Mark that the HDF5 object ``h5object`` (dataset or group) contains data with an
    iteration step (stride) of iter_step).'''
    h5group.attrs['iter_step'] = iter_step
    
def get_iter_step(h5group):
    '''Read back iteration step (stride) written by ``stamp_iter_step``'''
    return int(h5group.attrs['iter_step'])

def check_iter_range_least(h5object, iter_start, iter_stop):
    '''Return True if the iteration range [iter_start, iter_stop) is
    the same as or entirely contained within the iteration range stored
    on ``h5object``.'''
    obj_iter_start, obj_iter_stop = get_iter_range(h5object)
    return (obj_iter_start <= iter_start and obj_iter_stop >= iter_stop)

def check_iter_range_equal(h5object, iter_start, iter_stop):
    '''Return True if the iteration range [iter_start, iter_stop) is
    the same as the iteration range stored on ``h5object``.'''
    obj_iter_start, obj_iter_stop = get_iter_range(h5object)    
    return (obj_iter_start == iter_start and obj_iter_stop == iter_stop)

def get_iteration_entry(h5object, n_iter):
    '''Create a slice for data corresponding to iteration ``n_iter`` in ``h5object``.'''
    obj_iter_start, obj_iter_stop = get_iter_range(h5object)
    if n_iter < obj_iter_start or n_iter >= obj_iter_stop:
        raise IndexError('data for iteration {} not available in dataset {!r}'.format(n_iter, h5object))
    return numpy.index_exp[n_iter-obj_iter_start]

def get_iteration_slice(h5object, iter_start, iter_stop=None, iter_stride=None):
    '''Create a slice for data corresponding to iterations [iter_start,iter_stop),
    with stride iter_step, in the given ``h5object``.'''
    obj_iter_start, obj_iter_stop = get_iter_range(h5object)
    
    if iter_stop is None: iter_stop = iter_start+1
    if iter_stride is None: iter_stride = 1
    
    if iter_start < obj_iter_start:
        raise IndexError('data for iteration {} not available in dataset {!r}'.format(iter_start, h5object))
    elif iter_start > obj_iter_stop:
        raise IndexError('data for iteration {} not available in dataset {!r}'.format(iter_stop, h5object))
    
    start_index = iter_start - obj_iter_start
    stop_index = iter_stop - obj_iter_start
    return numpy.index_exp[start_index:stop_index:iter_stride]    



        
    
###
# Axis label metadata
###
def label_axes(h5object, labels, units=None):
    '''Stamp the given HDF5 object with axis labels. This stores the axis labels
    in an array of strings in an attribute called ``axis_labels`` on the given
    object. ``units`` if provided is a corresponding list of units.'''
    
    if len(labels) != len(h5object.shape):
        raise ValueError('number of axes and number of labels do not match')
    
    if units is None: units = []
    
    if len(units) and len(units) != len(labels):
        raise ValueError('number of units labels does not match number of axes')
    
    h5object.attrs['axis_labels'] = numpy.array(map(str,labels))
     
    if len(units):
        h5object.attrs['axis_units'] = numpy.array(map(str,units))

NotGiven = object()
def _get_one_attr(h5object, namelist, default=NotGiven):
    attrs = dict(h5object.attrs)
    for name in namelist:
        try:
            return attrs[name]
        except KeyError:
            pass
    else:
        if default is NotGiven:
            raise KeyError('no such key')
        else:
            return default
        
class WESTPAH5File(h5py.File):
    '''Generalized input/output for WESTPA simulation (or analysis) data.'''
    
    default_iter_prec = 8
    _this_fileformat_version = 8
        
    def __init__(self, *args, **kwargs):
        
        # These values are used for creating files or reading files where this
        # data is not stored. Otherwise, values stored as attributes on the root
        # group are used instead.
        arg_iter_prec = kwargs.pop('westpa_iter_prec', self.default_iter_prec)
        arg_fileformat_version = kwargs.pop('westpa_fileformat_version', self._this_fileformat_version)
        arg_creating_program = kwargs.pop('creating_program', None)
        
        # Initialize h5py file
        super(WESTPAH5File,self).__init__(*args, **kwargs)
        
        # Try to get iteration precision and I/O class version
        h5file_iter_prec = _get_one_attr(self, ['westpa_iter_prec', 'west_iter_prec', 'wemd_iter_prec'], None)
        h5file_fileformat_version = _get_one_attr(self,
                                                  ['westpa_fileformat_version', 
                                                   'west_file_format_version', 
                                                   'wemd_file_format_version'],
                                                  None)
        
        self.iter_prec = h5file_iter_prec if h5file_iter_prec is not None else arg_iter_prec
        self.fileformat_version = h5file_fileformat_version if h5file_fileformat_version is not None else arg_fileformat_version
        
        # Ensure that file format attributes are stored, if the file is writable
        if self.mode == 'r+':
            self.attrs['westpa_iter_prec'] = self.iter_prec
            self.attrs['westpa_fileformat_version'] = self.fileformat_version
            if arg_creating_program:
                stamp_creator_data(self, creating_program=arg_creating_program)
    
    # Iteration groups
    
    def iter_object_name(self, n_iter, prefix='', suffix=''):
        '''Return a properly-formatted per-iteration name for iteration
        ``n_iter``. (This is used in create/require/get_iter_group, but may 
        also be useful for naming datasets on a per-iteration basis.)'''
        return '{prefix}iter_{n_iter:0{prec}d}{suffix}'\
                .format(n_iter=n_iter,prefix=prefix,suffix=suffix,prec=self.iter_prec)
        
    def create_iter_group(self, n_iter, group=None):
        '''Create a per-iteration data storage group for iteration number ``n_iter``
        in the group ``group`` (which is '/iterations' by default).'''
        
        if group is None:
            group = self.require_group('/iterations')
        return group.create_group(self.iter_object_name(n_iter))
                
    def require_iter_group(self, n_iter, group=None):
        '''Ensure that a per-iteration data storage group for iteration number ``n_iter``
        is available in the group ``group`` (which is '/iterations' by default).'''
        if group is None:
            group = self.require_group('/iterations')
        return group.require_group(self.iter_object_name(n_iter))
        
    def get_iter_group(self, n_iter, group=None):
        '''Get the per-iteration data group for iteration number ``n_iter`` from within
        the group ``group`` ('/iterations' by default).'''
        if group is None:
            group = self['/iterations']
        return group[self.iter_object_name(n_iter)]
    
### Generalized WE dataset access classes
    
class DSSpec:
    '''Generalized WE dataset access'''
    
    def get_iter_data(self, n_iter, seg_slice=index_exp[:]):
        raise NotImplementedError
    
    def get_segment_data(self, n_iter, seg_id):
        return self.get_iter_data(n_iter)[seg_id]
    
    def __getstate__(self):
        d = dict(self.__dict__)
        if '_h5file' in d:
            d['_h5file'] = None
        return d
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        
    
class FileLinkedDSSpec(DSSpec):
    '''Provide facilities for accessing WESTPA HDF5 files, including auto-opening and the ability
    to pickle references to such files for transmission (through, e.g., the work manager), provided
    that the HDF5 file can be accessed by the same path on both the sender and receiver.'''
    def __init__(self, h5file_or_name):
        self._h5file = None
        self._h5filename = None
        
        try:
            self._h5filename = os.path.abspath(h5file_or_name.filename)
        except AttributeError:
            self._h5filename = h5file_or_name
            self._h5file = None
        else:
            self._h5file = h5file_or_name
            
    @property
    def h5file(self):
        if self._h5file is None:
            self._h5file = WESTPAH5File(self._h5filename, 'r')
        return self._h5file
    
class SingleDSSpec(FileLinkedDSSpec):
    @classmethod
    def from_string(cls, dsspec_string, default_h5file):
        alias = None
        
        h5file = default_h5file
        fields = dsspec_string.split(',')
        dsname = fields[0]
        
        for field in (field.strip() for field in fields[1:]):
            k,v = field.split('=')
            k = k.lower()
            if k == 'alias':
                alias = v
            elif k == 'slice':
                try:
                    slice = eval('numpy.index_exp' + v)
                except SyntaxError:
                    raise SyntaxError('invalid index expression {!r}'.format(v))
            elif k == 'file':
                h5file = WESTPAH5File(v, 'r')
            else:
                raise ValueError('invalid dataset option {!r}'.format(k))
            
        return cls(h5file, dsname, alias, slice)
        
    def __init__(self, h5file_or_name, dsname, alias=None, slice=None):
        FileLinkedDSSpec.__init__(self,h5file_or_name)
        self.dsname = dsname
        self.alias = alias or dsname
        self.slice = numpy.index_exp[slice] if slice else None
    
class SingleIterDSSpec(SingleDSSpec):
    def get_iter_data(self, n_iter, seg_slice=index_exp[:]):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[self.dsname][seg_slice + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[self.dsname][seg_slice]
    
class SingleSegmentDSSpec(SingleDSSpec):
    def get_iter_data(self, n_iter, seg_slice=index_exp[:]):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[self.dsname][seg_slice + index_exp[:] + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[self.dsname][seg_slice + index_exp[:,:]]
        
    def get_segment_data(self, n_iter, seg_id):
        if self.slice:
            return self.h5file.get_iter_group(n_iter)[numpy.index_exp[seg_id,:] + self.slice]
        else:
            return self.h5file.get_iter_group(n_iter)[seg_id,:]
        
class FnDSSpec(FileLinkedDSSpec):
    def __init__(self, h5file_or_name, fn):
        FileLinkedDSSpec.__init__(self,h5file_or_name)
        self.fn = fn
        
    def get_iter_data(self, n_iter, seg_slice=index_exp[:]):
        return self.fn(n_iter, self.h5file.get_iter_group(n_iter))[seg_slice]

        
class MultiDSSpec(DSSpec):
    def __init__(self, dsspecs):
        self.dsspecs = dsspecs
    
    def get_iter_data(self, n_iter, seg_slice=index_exp[:]):
        datasets = [dsspec.get_iter_data(n_iter) for dsspec in self.dsspecs]
          
        ncols = 0 
        nsegs = None
        npts = None
        for iset, dset in enumerate(datasets):
            if nsegs is None:
                nsegs = dset.shape[0]
            elif dset.shape[0] != nsegs:
                raise TypeError('dataset {} has incorrect first dimension (number of segments)'.format(self.dsspecs[iset]))
            if npts is None:
                npts = dset.shape[1]
            elif dset.shape[1] != npts:
                raise TypeError('dataset {} has incorrect second dimension (number of time points)'.format(self.dsspecs[iset]))
            
            if dset.ndim < 2:
                # scalar per segment or scalar per iteration
                raise TypeError('dataset {} has too few dimensions'.format(self.dsspecs[iset]))
            elif dset.ndim > 3:
                # array per timepoint
                raise TypeError('dataset {} has too many dimensions'.format(self.dsspecs[iset]))
            elif dset.ndim == 2:
                # scalar per timepoint
                ncols += 1
            else:
                # vector per timepoint
                ncols += dset.shape[-1]
        
        output_dtype = numpy.result_type(*[ds.dtype for ds in datasets])
        output_array = numpy.empty((nsegs, npts, ncols), dtype=output_dtype)
        
        ocol = 0
        for iset, dset in enumerate(datasets):
            if dset.ndim == 2:
                output_array[:,:,ocol] = dset[...]
                ocol += 1
            elif dset.ndim == 3:
                output_array[:,:,ocol:(ocol+dset.shape[-1])] = dset[...]
                ocol += dset.shape[-1]
        
        return output_array[seg_slice]
            
        
        
        
        
        
        
        
        
        