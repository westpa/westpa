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


console_scripts_core = [
    'w_fork = westpa.cli.core.w_fork:entry_point',
    'w_states = westpa.cli.core.w_states:entry_point',
    'w_run = westpa.cli.core.w_run:entry_point',
    'w_truncate = westpa.cli.core.w_truncate:entry_point',
    'w_init = westpa.cli.core.w_init:entry_point',
    'w_succ = westpa.cli.core.w_succ:entry_point',
]


console_scripts_tools = [
    'w_direct = westpa.cli.tools.w_direct:entry_point',
    'w_dumpsegs = westpa.cli.tools.w_dumpsegs:entry_point',
    'w_postanalysis_matrix = westpa.cli.tools.w_postanalysis_matrix:entry_point',
    'w_select = westpa.cli.tools.w_select:entry_point',
    'w_ntop = westpa.cli.tools.w_ntop:entry_point',
    'w_kinavg = westpa.cli.tools.w_kinavg:entry_point',
    'w_eddist = westpa.cli.tools.w_eddist:entry_point',
    'w_assign = westpa.cli.tools.w_assign:entry_point',
    'w_trace = westpa.cli.tools.w_trace:entry_point',
    'w_crawl = westpa.cli.tools.w_crawl:entry_point',
    'w_kinetics = westpa.cli.tools.w_kinetics:entry_point',
    'w_fluxanl = westpa.cli.tools.w_fluxanl:entry_point',
    'w_reweight = westpa.cli.tools.w_reweight:entry_point',
    'w_pdist = westpa.cli.tools.w_pdist:entry_point',
    'w_ipa = westpa.cli.tools.w_ipa:entry_point',
    'w_bins = westpa.cli.tools.w_bins:entry_point',
    'w_stateprobs = westpa.cli.tools.w_stateprobs:entry_point',
    'w_postanalysis_reweight = westpa.cli.tools.w_postanalysis_reweight:entry_point',
    'ploterr = westpa.cli.tools.ploterr:entry_point',
    'plothist = westpa.cli.tools.plothist:entry_point',
    'w_multi_west = westpa.cli.tools.w_multi_west:entry_point',
    'w_red = westpa.cli.tools.w_red:entry_point',
    'w_timings = westpa.cli.tools.w_timings:entry_point',
]


console_scripts = console_scripts_core + console_scripts_tools


metadata = dict(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    entry_points={'console_scripts': console_scripts},
)


if __name__ == '__main__':
    metadata['ext_modules'] = extensions()
    setup(**metadata)
