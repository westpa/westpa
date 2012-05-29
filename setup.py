from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
numpy_include = numpy.get_include()

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("fasthist._fasthist", 
                               ["fasthist/_fasthist.pyx"], #,"fasthist/_fasthist_supp.c"], 
                               include_dirs=['.', numpy_include],
                               extra_compile_args=['-O2'])])

