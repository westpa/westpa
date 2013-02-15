'''Miscellaneous routines to help with HDF5 input and output of WEST-related data.'''

import sys, getpass, socket, time
import numpy, h5py

#
# Helper functions
# 
def calc_chunksize(shape, dtype, max_chunksize=262144):
    '''Calculate a chunk size for HDF5 data, anticipating that access will slice
    along lower dimensions sooner than higher dimensions.'''
        
    chunk_shape = list(shape)
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

def stamp_iter_range(h5object, start_iter, stop_iter):
    '''Mark that the HDF5 object ``h5object`` (dataset or group) contains data from iterations
    start_iter <= n_iter < stop_iter.'''
    attrs = h5object.attrs
    attrs['iter_start'] = start_iter
    attrs['iter_stop']  = stop_iter
    
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
    