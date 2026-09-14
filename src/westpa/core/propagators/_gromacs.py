import os
import shutil
import subprocess

import numpy as np

from .base import SerialPropagator
from ..state import State


class GROMACSPropagator(SerialPropagator):
    """`GROMACS <https://www.gromacs.org/>`_ molecular dynamics propagator.

    To create an initial state for this propagator, pass the absolute path of
    a GROMACS coordinate file to the :class:`State` constructor's `file` parameter::

        state = westpa.State(file=os.path.abspath('conf.gro'))

    Parameters
    ----------
    top_file : path-like
        Topology file.
    mdp_file : path-like
        MD parameter file. The ``ld-seed`` option is dynamically overridden
        for each input segment.
    grompp_args : sequence of str, optional
        Optional arguments to pass to ``gmx grompp``. The ``-p``, ``-c`` and
        ``-o`` arguments are reserved and cannot be overridden.
    mdrun_args : sequence of str, optional
        Optional arguments to pass to ``gmx mdrun``. The ``-s`` argument is
        reserved and cannot be overridden.
    **kwargs
        Keyword arguments to pass to the :class:`SerialPropagator` base class.

    Examples
    --------
    >>> import westpa
    >>> propagator = westpa.GROMACSPropagator('topol.top', 'grompp.mdp')

    """

    def __init__(
        self,
        top_file,
        mdp_file,
        grompp_args=None,
        mdrun_args=None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        if not shutil.which('gmx'):
            raise RuntimeError("couldn't find 'gmx' executable")

        grompp_args = list(grompp_args) if grompp_args is not None else []
        if {'-p', '-c', '-o'}.intersection(grompp_args):
            raise ValueError("invalid 'grompp_args' input: the -p, -c, and -o flags are not supported")

        mdrun_args = list(mdrun_args) if mdrun_args is not None else []
        if '-s' in mdrun_args:
            raise ValueError("invalid 'mdrun_args' input: the -s flag is not supported")

        self.top_file = os.path.abspath(top_file)
        self.mdp_file = os.path.abspath(mdp_file)
        self.grompp_args = grompp_args
        self.mdrun_args = mdrun_args

    def propagate(self, segment, rng):
        segment_dir = self.make_segment_dir(segment)

        # copy .mdp file to segment directory
        try:
            idx = self.grompp_args.index('-f')
        except ValueError:
            mdp_file = os.path.join(segment_dir, 'grompp.mdp')  # gmx default
        else:
            mdp_file = os.path.join(segment_dir, self.args[idx + 1])
        shutil.copyfile(self.mdp_file, mdp_file)

        # override 'ld-seed' option
        ld_seed = rng.integers(2**16, dtype=np.uint16)
        with open(mdp_file, 'a') as f:
            f.write(f'\nld-seed = {ld_seed}\n')

        # assemble grompp and mdrun commands
        grompp_cmd = ['gmx', 'grompp', '-p', self.top_file, '-c', segment.initial_state.file, *self.grompp_args]
        mdrun_cmd = ['gmx', 'mdrun', *self.mdrun_args]

        # call grompp and mdrun
        for cmd in grompp_cmd, mdrun_cmd:
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
            idx = self.mdrun_args.index('-c')
        except ValueError:
            filename = 'confout.gro'  # gmx default
        else:
            filename = self.mdrun_args[idx + 1]
        segment.final_state = State(file=os.path.join(segment_dir, filename))

        return segment
