'''Miscellaneous routines to help with HDF5 input and output of WEST-related data.'''

import sys, getpass, socket, time
import numpy, h5py

def create_hdf5_group(parent_group, groupname, replace=False, creating_program=None):
    '''Create (or delete and recreate) and HDF5 group named ``groupname`` within
    the enclosing Group (object) ``parent_group``. If ``replace`` is True, then
    the group is replaced if present; if False, then an error is raised if the
    group is present. After the group is created, HDF5 attributes are set recording
    the following:
    
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
    
    if replace:
        try:
            del parent_group[groupname]
        except KeyError:
            pass
        
    newgroup = parent_group.create_group(groupname)
    now = time.time()
    attrs = newgroup.attrs
    
    attrs['creation_program'] = creating_program or sys.argv[0] or 'unknown program'
    attrs['creation_user'] = getpass.getuser()
    attrs['creation_hostname'] = socket.gethostname()
    attrs['creation_unix_time'] = now
    attrs['creation_time'] = time.strftime('%c', time.localtime(now))
    
    return newgroup 

def label_axes(h5object, labels, units=None):
    '''Stamp the given HDF5 object with axis labels. This stores the axis labels
    in an array of strings in an attribute called ``axis_labels`` on the given
    object.'''
    
    if len(labels) != len(h5object.shape):
        raise ValueError('number of axes and number of labels do not match')
    
    if units is None: units = []
    
    if len(units) and len(units) != len(labels):
        raise ValueError('number of units labels does not match number of axes')
    
    h5object.attrs['axis_labels'] = numpy.array(map(str,labels))
     
    if len(units):
        h5object.attrs['axis_units'] = numpy.array(map(str,units))
    