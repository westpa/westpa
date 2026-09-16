import logging
import os
import shutil
import subprocess

import f90nml
import numpy as np

from .base import SerialPropagator
from ..state import State

logger = logging.getLogger(__name__)


class AmberPropagator(SerialPropagator):
    """`Amber <https://ambermd.org/>`_ molecular dynamics propagator.

    To create an initial state for this propagator, pass the absolute path of
    an Amber coordinate file to the :class:`State` constructor's `file` parameter::

        state = westpa.State(file=os.path.abspath('inpcrd'))

    Parameters
    ----------
    prmtop_file : path-like
        Parameter-topology file.
    mdin_file : path-like
        MD input file. The ``ig``, ``irest``, and ``ntx`` options are
        dynamically overridden for each input segment.
    engine : str, optional
        Amber program to execute. Defaults to ``'pmemd.cuda'``, ``'pmemd'``,
        or ``'sander'``, whichever is found first.
    args : sequence of str, optional
        Optional arguments to pass to `engine`.
        The ``-p`` and ``-c`` arguments are reserved and cannot be overridden.
    **kwargs
        Keyword arguments to pass to the :class:`SerialPropagator` base class.

    Examples
    --------

    >>> import westpa
    >>> propagator = westpa.AmberPropagator('prmtop', 'mdin')

    """

    def __init__(
        self,
        prmtop_file,
        mdin_file,
        engine=None,
        args=None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        if engine is None:
            for engine in ('pmemd.cuda', 'pmemd', 'sander'):
                if shutil.which(engine):
                    break
            else:
                raise RuntimeError('default Amber engine (pmemd.cuda, pmemd, or sander) not found')
        elif not shutil.which(engine):
            raise RuntimeError(f"couldn't find {engine!r} executable")

        args = list(args) if args is not None else []
        if '-p' in args or '-c' in args:
            raise ValueError("invalid 'args' input: the -p and -c flags are not supported")

        self.prmtop_file = os.path.abspath(prmtop_file)
        self.mdin_file = os.path.abspath(mdin_file)
        self.engine = engine
        self.args = args

    def propagate(self, segment, rng):
        segment_dir = self.make_segment_dir(segment)

        # read user-provided mdin file
        with open(self.mdin_file, 'r') as f:
            title_line = f.readline()
            mdin = f90nml.read(f)

        # override 'ig', 'irest', and 'ntx' options
        ig = rng.integers(2**16, dtype=np.uint16).item()
        if segment.initpoint_type == segment.InitPoint.NEWTRAJ:
            irest, ntx = 0, 1
        else:
            irest, ntx = 1, 5
        mdin['cntrl'].update(ig=ig, irest=irest, ntx=ntx)

        # write mdin file to segment directory
        try:
            idx = self.args.index('-i')
        except ValueError:
            mdin_file = os.path.join(segment_dir, 'mdin')  # Amber default
        else:
            mdin_file = os.path.join(segment_dir, self.args[idx + 1])
        with open(mdin_file, 'w') as f:
            f.write(title_line)
            f90nml.write(mdin, f)

        # run amber
        cmd = [self.engine, '-p', self.prmtop_file, '-c', segment.initial_state.file, *self.args]
        completed_process = subprocess.run(
            cmd,
            cwd=segment_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        with open(os.path.join(segment_dir, 'stderr.txt'), 'a') as f:
            f.write(completed_process.stderr)

        # retrieve the final state
        try:
            idx = self.args.index('-r')
        except ValueError:
            filename = 'restrt'  # Amber default
        else:
            filename = self.args[idx + 1]
        segment.final_state = State(file=os.path.join(segment_dir, filename))

        return segment
