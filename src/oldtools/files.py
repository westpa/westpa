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
