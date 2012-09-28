from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
numpy_include = numpy.get_include()

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("west.binning._assign", 
                               ["west/binning/_assign.pyx"], 
                               include_dirs=['.', numpy_include],                               
                               # hack-ish; included since my dev box has trouble
                               extra_compile_args=['-O3']),
                     Extension("west.kinetics._kinetics", 
                                ["west/kinetics/_kinetics.pyx"], 
                                include_dirs=['.', numpy_include],
                                extra_compile_args=['-O3']),
                     ])
