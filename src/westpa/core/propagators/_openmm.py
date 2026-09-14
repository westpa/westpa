import copy
import logging
import os
from dataclasses import dataclass

import openmm.app

from .base import SerialPropagator
from westpa.core.state import State

logger = logging.getLogger(__name__)


@dataclass
class Report:
    reporter_type: type
    filename: str
    report_interval: int
    options: dict


class OpenMMPropagator(SerialPropagator):
    """`OpenMM <https://openmm.org/>`_ molecular dynamics propagator.

    To create an initial state for this propagator,
    `save <https://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html#openmm.app.simulation.Simulation.saveState>`_
    an OpenMM state to an XML file, and pass the absolute file path to the
    :class:`State` constructor's `file` parameter::

        state = westpa.State(file=os.path.abspath('state.xml'))

    Parameters
    ----------
    topology : openmm.app.Topology
        Molecular topology (chains, residues, atoms, and bonds).
    system : openmm.System
        OpenMM System object created by applying a force field to `topology`.
    integrator : openmm.Integrator
        Integrator to use for simulating the system.
    steps : int
        Number of time steps to simulate for each segment.
    platform : openmm.Platform, optional
        Platform to use for calculations.
    platform_properties : Mapping[str, str], optional
        Platform-specific properties to pass to the simulation context.
    final_state_filename : str, optional
        Name of the XML file used to store the final state of a segment.
        Defaults to ``'final_state.xml'``.
    **kwargs
        Keyword arguments to pass to the :class:`SerialPropagator` base class.

    Examples
    --------

    Create a propagator that runs 2 ps of Langevin dynamics at 300 K:

    >>> import westpa
    >>> import openmm.app
    >>> from openmm import unit
    >>> pdb = openmm.app.PDBFile('topology.pdb')
    >>> forcefield = openmm.app.ForceField('amber14-all.xml')
    >>> propagator = westpa.OpenMMPropagator(
    ...     topology=pdb.topology,
    ...     system=forcefield.createSystem(
    ...         pdb.topology,
    ...         nonbondedMethod=openmm.app.PME,
    ...         constraints=openmm.app.HBonds,
    ...     ),
    ...     integrator=openmm.LangevinIntegrator(
    ...         temperature=300 * unit.kelvin,
    ...         frictionCoeff=1 / unit.picosecond,
    ...         stepSize=2 * unit.femtosecond,
    ...     ),
    ...     steps=1000,
    ... )

    Add a reporter that writes the positions of the first 22 atoms every 100 steps:

    >>> propagator.add_reporter(
    ...     openmm.app.XTCReporter,
    ...     filename='traj.xtc',
    ...     report_interval=100,
    ...     options={'atomSubset': list(range(22))},
    ... )

    Add a reporter that writes the kinetic and potential energy every 500 steps:

    >>> propagator.add_reporter(
    ...     openmm.app.StateDataReporter,
    ...     filename='log.csv',
    ...     report_interval=500,
    ...     options={'kineticEnergy': True, 'potentialEnergy': True},
    ... )

    Notes
    -----

    To use this propagator, the `OpenMM <https://pypi.org/project/OpenMM/>`_
    package must be installed, for example:

    .. code-block:: shell

       conda install conda-forge::openmm

    or (since OpenMM 8.1.1):

    .. code-block:: shell

       pip install openmm

    """

    DEFAULT_FINAL_STATE_FILENAME = 'final_state.xml'

    def __init__(
        self,
        topology,
        system,
        integrator,
        steps,
        platform=None,
        platform_properties=None,
        final_state_filename=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.topology = topology
        self.system = system
        self.integrator = integrator
        self.steps = steps
        self.platform = platform
        self.platform_properties = platform_properties
        self.final_state_filename = final_state_filename or self.DEFAULT_FINAL_STATE_FILENAME
        self._reports = []

    def add_reporter(self, reporter_type, filename, report_interval, options=None):
        """Add a reporter to save time series data for each segment.

        Parameters
        ----------
        reporter_type : type
            Class that implements the OpenMM reporter protocol. It must be
            possible to create a reporter instance by calling
            ``reporter_type(filename, report_interval, **options)``.
        filename : str
            Name of the file to write output to.
        report_interval : int
            Interval (in time steps) at which to write frames.
        options : Mapping[str, Any], optional
            Optional keyword arguments to pass to the `reporter_type` constructor.

        """
        report = Report(reporter_type, filename, report_interval, options or {})
        self._reports.append(report)

    def propagate(self, segment, rng):
        integrator = copy.copy(self.integrator)  # integrators can only bind one context

        if hasattr(integrator, 'setRandomNumberSeed'):
            integrator.setRandomNumberSeed(rng.integers(low=1, high=2**31))
        for force in self.system.getForces():
            if hasattr(force, 'setRandomNumberSeed'):
                force.setRandomNumberSeed(rng.integers(low=1, high=2**31))

        simulation = openmm.app.Simulation(
            self.topology,
            self.system,
            integrator,
            platform=self.platform,
            platformProperties=self.platform_properties,
            state=segment.initial_state.file,
        )

        segment_dir = self.make_segment_dir(segment)

        for report in self._reports:
            file = os.path.join(segment_dir, report.filename)
            reporter = report.reporter_type(file, report.report_interval, **report.options)
            simulation.reporters.append(reporter)

        simulation.step(self.steps)

        final_state_file = os.path.join(segment_dir, self.final_state_filename)
        simulation.saveState(final_state_file)
        segment.final_state = State(file=final_state_file)

        return segment
