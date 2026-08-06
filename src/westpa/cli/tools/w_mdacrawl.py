import importlib
import inspect
import logging
import sys

import numpy as np

from westpa.tools import WESTTool
import MDAnalysis as mda
from westpa.core import mdanalysis_core  # noqa
from westpa.core.mdanalysis_core import save_to_west_h5

log = logging.getLogger('w_mdacrawl')


class WMDACrawl(WESTTool):
    prog = 'w_mdacrawl'
    description = '''\
Run an MDAnalysis analysis on a WESTPA west.h5 file.
The analysis is run on the Universe created by WESTPAParser and WESTPAReader.
Results can optionally be saved back into the HDF5 as auxdata.

Examples:
  w_mdacrawl west.h5 MDAnalysis.analysis.rms.RMSD --save rmsd --select "name CA" --column 2
  w_mdacrawl west.h5 MDAnalysis.analysis.rms.RMSF --save rmsf --select "protein"
  w_mdacrawl west.h5 MDAnalysis.analysis.rms.RMSD --save rmsd -j 4
-------------------------------------------------------------------------------------------
'''

    def __init__(self):
        super().__init__()
        self.analysis_cls = None
        self.universe = None
        self.results = None
        self.args = None

    def add_args(self, parser):
        parser.add_argument(
            'west_h5',
            help='Path to the WESTPA west.h5 file',
        )
        parser.add_argument(
            'analysis',
            help=(
                'Dotted path to the MDAnalysis analysis class '
                '(e.g., MDAnalysis.analysis.rms.RMSD) or a user-defined '
                'module.ClassName'
            ),
        )
        parser.add_argument(
            '--column',
            type=int,
            default=None,
            help='Extract a specific column index from 2D results (e.g., 2 for RMSD values)',
        )
        parser.add_argument(
            '--save',
            metavar='DATASET_NAME',
            default=None,
            help=('Save the results back into west.h5 under ' 'iterations/iter_XXXXXXXX/auxdata/DATASET_NAME.'),
        )
        parser.add_argument(
            '--select',
            metavar='SELECTION',
            default=None,
            help='MDAnalysis atom selection string (e.g., "name CA", "protein")',
        )
        parser.add_argument(
            '-j',
            '--n-workers',
            type=int,
            default=1,
            metavar='N',
            help='Number of parallel workers for multiprocessing backend (default: 1 = serial)',
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            default=False,
            help='Overwrite the dataset if it already exists in the auxdata dataset',
        )

    def process_args(self, args):
        super().process_args(args)
        self.args = args

        log.info(f'Importing analysis class: {self.args.analysis}')
        self.analysis_cls = self._import_analysis_class(self.args.analysis)
        log.info(f'Loaded: {self.analysis_cls}')

    def go(self):
        log.info(f'Loading Universe from {self.args.west_h5}')
        self.universe = mda.Universe(self.args.west_h5, format='WESTPA')
        log.info(f'Universe loaded: {self.universe.atoms.n_atoms} atoms, {len(self.universe.trajectory)} frames')

        analysis = self._run_analysis(self.analysis_cls, self.universe, self.args.select, self.args.n_workers, {})
        log.info('Analysis complete')

        # Extract and slice results
        self.results = self._extract_results(analysis)
        if self.args.column is not None:
            if self.results.ndim > 1:
                log.info(f'Slicing column {self.args.column} from results.')
                self.results = self.results[:, self.args.column]
            else:
                log.warning('--column ignored: Results are already a 1D array.')

        log.info(f'Results shape: {self.results.shape}, dtype: {self.results.dtype}')
        print(f'\n{"="*60}')
        print(f'Analysis : {self.analysis_cls.__name__}')
        print(f'Frames   : {len(self.universe.trajectory)}')
        print(f'Results  : shape={self.results.shape}, dtype={self.results.dtype}')

        if self.results.ndim <= 2 and self.results.shape[0] <= 10:
            print(f'Values   :\n{self.results}')
        else:
            print(f'First 5  :\n{self.results[:5]}')
            print(f'Last 5   :\n{self.results[-5:]}')
        print(f'{"="*60}\n')

        # Release the HDF5 read-lock before triggering any potential save events
        self.universe.trajectory.close()

        if self.args.save:
            log.info(f'Saving results as auxdata/{self.args.save} (overwrite={self.args.overwrite})')
            save_to_west_h5(
                self.universe,
                self.results,
                self.args.save,
                west_h5_path=self.args.west_h5,
                overwrite=self.args.overwrite,
            )
            log.info(f'Successfully saved to {self.args.west_h5}::auxdata/{self.args.save}')
            print(f'Results saved to {self.args.west_h5} under auxdata/{self.args.save}')
        else:
            print(f'Results were NOT saved (use --save DATASET_NAME to write to {self.args.west_h5})')

    # Helper Methods
    def _import_analysis_class(self, dotted_path):
        parts = dotted_path.rsplit('.', 1)
        if len(parts) != 2:
            raise ValueError(f"analysis must be specified as 'module.ClassName', got {dotted_path!r}")

        module_path, class_name = parts

        if '.' not in module_path:
            sys.path.insert(0, '.')

        try:
            module = importlib.import_module(module_path)
        except ModuleNotFoundError as e:
            raise ImportError(
                f"Could not import module {module_path!r}. "
                f"Make sure the package is installed or the script is in the current directory."
            ) from e

        try:
            cls = getattr(module, class_name)
        except AttributeError:
            raise AttributeError(
                f"Module {module_path!r} has no attribute {class_name!r}. "
                f"Available: {[x for x in dir(module) if not x.startswith('_')]}"
            )
        return cls

    def _run_analysis(self, analysis_cls, universe, select, n_workers, extra_kwargs):
        sig = inspect.signature(analysis_cls.__init__)
        params = [p for p in list(sig.parameters.keys()) if p != 'self']
        init_kwargs = {}

        if select:
            if 'select' in params:
                init_kwargs['select'] = select
            elif 'atomgroup' in params or (params and params[0] in ('atomgroup', 'ag')):
                init_kwargs[params[0]] = universe.select_atoms(select)
            else:
                init_kwargs[params[0]] = universe.select_atoms(select)

        init_kwargs.update(extra_kwargs)
        needs_reference = any(p in params for p in ('reference', 'ref'))

        if needs_reference:
            if select:
                analysis = analysis_cls(universe, universe, select=select, **extra_kwargs)
            else:
                analysis = analysis_cls(universe, universe, **extra_kwargs)
        elif select:
            if 'select' in params:
                analysis = analysis_cls(universe, select=select, **extra_kwargs)
            else:
                ag = universe.select_atoms(select)
                analysis = analysis_cls(ag, **extra_kwargs)
        else:
            analysis = analysis_cls(universe, **extra_kwargs)

        if n_workers and n_workers > 1:
            log.info(f'Running with multiprocessing backend, n_workers={n_workers}')
            analysis.run(backend='multiprocessing', n_workers=n_workers)
        else:
            log.info('Running with serial backend')
            analysis.run(backend='serial')

        return analysis

    def _extract_results(self, analysis):
        results = analysis.results
        # Try common result attribute names in order of preference
        for attr in ['rmsd', 'rmsf', 'rdf', 'timeseries', 'results']:
            if hasattr(results, attr):
                val = getattr(results, attr)
                if isinstance(val, np.ndarray):
                    return val

        # Results is often a dict-like object (MDAnalysis Results class)
        if hasattr(results, 'items'):
            for key, val in results.items():
                if isinstance(val, np.ndarray) and val.size > 0:
                    log.info(f'Auto-detected result attribute: {key}')
                    return val

        # Fallback to __dict__ just in case
        if hasattr(results, '__dict__'):
            for key, val in results.__dict__.items():
                if isinstance(val, np.ndarray) and val.size > 0:
                    log.info(f'Auto-detected result attribute: {key}')
                    return val

        raise RuntimeError(
            f"Could not auto-detect results from {type(analysis).__name__}. "
            f"Available result attributes/keys: {list(getattr(results, 'keys', lambda: vars(results).keys())())}"
        )


def entry_point():
    WMDACrawl().main()


if __name__ == '__main__':
    entry_point()
