Python API
==========

The public API is exported at the package level::

   import westpa


Data types
----------

.. autosummary::
   :nosignatures:

   ~westpa.State
   ~westpa.Segment
   ~westpa.Bin
   ~westpa.Source
   ~westpa.Sink


Simulations
-----------

.. autosummary::
   :nosignatures:

   ~westpa.Simulation


Analysis
--------

.. autosummary::
   :nosignatures:

   ~westpa.TrajectoryTree
   ~westpa.TrajectoryTreeView


Propagators
-----------

Protocol
~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.Propagator

Base classes
~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.SerialPropagator
   ~westpa.VectorizedPropagator

Built-in implementations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.AmberPropagator
   ~westpa.GROMACSPropagator
   ~westpa.OpenMMPropagator


Bin mappers
-----------

Protocol
~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.BinMapper

Built-in implementations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.RectilinearBinMapper
   ~westpa.MABBinMapper
   ~westpa.VoronoiBinMapper
   ~westpa.AdaptiveVoronoiBinMapper


Resamplers
----------

Protocol
~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.Resampler

Base class
~~~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.ResamplerBase

Built-in implementations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   ~westpa.HuberKimResampler
   ~westpa.MultinomialResampler
   ~westpa.ResidualResampler
   ~westpa.StratifiedResampler
   ~westpa.SystematicResampler

Work managers
-------------

.. autosummary::
   :nosignatures:

   ~westpa.SerialWorkManager
   ~westpa.ThreadsWorkManager
   ~westpa.ProcessWorkManager
