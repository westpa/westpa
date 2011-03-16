from __future__ import division, print_function; __metaclass__ = type

import numpy

def load_npy_or_text(filename, **kwargs):
    '''Load an array from an existing .npy file, or read a text file and
    convert to a NumPy array.  In either case, return a NumPy array.  If a 
    pickled NumPy dataset is found, memory-map it read-only.  If the specified
    file does not contain a pickled NumPy array, attempt to read the file using
    numpy.loadtxt(filename, **kwargs).'''
    
    try:
        return numpy.load(filename, 'r')
    except IOError as e:
        if 'Failed to interpret' in str(e):
            pass
        else:
            raise
    
    return numpy.loadtxt(filename, **kwargs)