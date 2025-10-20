import sys

from setuptools import setup, Extension

import versioneer


def extensions():
    from Cython.Build import cythonize
    from numpy import get_include as np_get_include

    np_inc = np_get_include()

    common_cflags = [
        '-O3',
    ]

    fasthist_module = Extension(
        'westpa.fasthist._fasthist',
        sources=['src/westpa/fasthist/_fasthist.pyx'],
        include_dirs=[np_inc],
        extra_compile_args=common_cflags,
    )

    trajtree_module = Extension(
        'westpa.trajtree._trajtree',
        sources=['src/westpa/trajtree/_trajtree.pyx'],
        include_dirs=[np_inc],
        extra_compile_args=common_cflags,
    )

    mclib_module = Extension(
        'westpa.mclib._mclib', sources=['src/westpa/mclib/_mclib.pyx'], include_dirs=[np_inc], extra_compile_args=common_cflags
    )

    binning_module = Extension(
        'westpa.core.binning._assign',
        sources=['src/westpa/core/binning/_assign.pyx'],
        include_dirs=[np_inc],
        extra_compile_args=common_cflags,
    )

    kinetics_module = Extension(
        'westpa.core.kinetics._kinetics',
        sources=['src/westpa/core/kinetics/_kinetics.pyx'],
        include_dirs=[np_inc],
        extra_compile_args=common_cflags,
    )

    reweight_module = Extension(
        'westpa.core.reweight._reweight',
        sources=['src/westpa/core/reweight/_reweight.pyx'],
        include_dirs=[np_inc],
        extra_compile_args=common_cflags,
    )

    exts = [fasthist_module, trajtree_module, mclib_module, binning_module, kinetics_module, reweight_module]

    exts = cythonize(exts, language_level=sys.version_info[0])

    return exts


metadata = dict(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)


if __name__ == '__main__':
    metadata['ext_modules'] = extensions()
    setup(**metadata)
