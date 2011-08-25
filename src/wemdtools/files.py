from __future__ import division, print_function; __metaclass__ = type

import os, sys, warnings
import numpy
import wemd
from wemd.util.config_dict import ConfigDict

def load_npy_or_text(filename):
    '''Load an array from an existing .npy file, or read a text file and
    convert to a NumPy array.  In either case, return a NumPy array.  If a 
    pickled NumPy dataset is found, memory-map it read-only.  If the specified
    file does not contain a pickled NumPy array, attempt to read the file using
    numpy.loadtxt(filename, **kwargs).'''
    
    f = open(filename, 'rb')
    try:
        f.seek(0)
    except IOError:
        # Not seekable - assume a text stream
        return numpy.loadtxt(filename)
    else:
        f.close()
        
    
    # File is seekable
    try:
        # try to mmap it
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return numpy.load(filename, 'r')
        
    except IOError as e:
        if 'Failed to interpret' in str(e):
            pass
        else:
            raise
    
    return numpy.loadtxt(filename)

def sim_manager_from_args(args, filename_attr='datafile', status_stream=None, script_name=None):
    '''Construct and return a SimManager object from a named file or from the file named in 
    ``wemd.cfg``, according to command line arguments.  Also initializes the logging subsystem.'''
    
    if script_name is None:
        script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    
    datafile = getattr(args,filename_attr)    
    
    if datafile:
        # Load from the given datafile
        runtime_config = ConfigDict()
        runtime_config['data.h5file'] = datafile
        if status_stream:
            status_stream.write('Reading WEMD data from {}\n'.format(datafile))
    else:
        # Load from the file specified in wemd.cfg
        if status_stream:
            status_stream.write('Reading WEMD data from HDF5 file specified in {}\n'.format(args.run_config_file or 'wemd.cfg'))
        runtime_config = wemd.rc.read_config(args.run_config_file) 
        
    runtime_config.update_from_object(args)
    wemd.rc.config_logging(args, script_name)
    return wemd.rc.load_sim_manager(runtime_config)
