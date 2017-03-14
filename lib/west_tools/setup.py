# Copyright (C) 2017 Matthew C. Zwier, Joshua L. Adelman, and Lillian T. Chong.
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

