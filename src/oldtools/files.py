


import warnings
import numpy

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
