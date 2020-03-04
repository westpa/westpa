
from distutils.core import setup
from distutils.extension import Extension

import numpy
numpy_include = numpy.get_include()

try:
    from Cython.Distutils import build_ext
    use_cython = True
    cmdclass = {'build_ext': build_ext}
except ImportError:
    use_cython = False
    cmdclass = {}
finally:
    print('Regenerating C extension code using Cython: {}'.format(use_cython))

suffix = 'pyx' if use_cython else 'c'

setup(cmdclass = cmdclass,
      ext_modules = [Extension("fasthist._fasthist", 
                               ["fasthist/_fasthist.{}".format(suffix)], #,"fasthist/_fasthist_supp.c"], 
                               include_dirs=['.', numpy_include],
                               
                               # hack-ish; included since my dev box has trouble
                               extra_compile_args=['-O3']),
                     Extension("trajtree._trajtree",
                               ["trajtree/_trajtree.{}".format(suffix)],
                               include_dirs=['.', numpy_include],
                               extra_compile_args=['-O3']),
                     Extension("mclib._mclib",
                               ["mclib/_mclib.{}".format(suffix)],
                               include_dirs=['.', numpy_include],
                               extra_compile_args=['-O3']),
                     Extension("westpa.binning._assign",
                               ["westpa/binning/_assign.{}".format(suffix)],
                               include_dirs=['.', numpy_include],
                               # hack-ish; included since my dev box has trouble
                               extra_compile_args=['-O3']),
                     Extension("westpa.kinetics._kinetics",
                                ["westpa/kinetics/_kinetics.{}".format(suffix)],
                                include_dirs=['.', numpy_include],
                                extra_compile_args=['-O3']),
                     Extension("westpa.reweight._reweight",
                                ["westpa/reweight/_reweight.{}".format(suffix)],
                                include_dirs=['.', numpy_include],
                                extra_compile_args=['-O3']),])

