Python API
==========

Representation
--------------

.. autosummary::
   :caption: Representation
   :toctree:
   :nosignatures:

   ~westpa.State
   ~westpa.Segment
   ~westpa.Bin
   ~westpa.Source
   ~westpa.Sink


Simulation
----------

.. autosummary::
   :caption: Simulation
   :toctree:
   :nosignatures:

   ~westpa.Simulation

Dynamics propagation
--------------------

.. autosummary::
   :caption: Dynamics propagation
   :toctree:
   :nosignatures:

   ~westpa.Propagator
   ~westpa.SerialPropagator
   ~westpa.VectorizedPropagator
   ~westpa.AmberPropagator
   ~westpa.GROMACSPropagator
   ~westpa.OpenMMPropagator

Calculating progress coordinates
--------------------------------

.. autosummary::
   :caption: Calculating progress coordinates
   :toctree:
   :nosignatures:

   ~westpa.PCoordCalculator


Binning
-------

.. autosummary::
   :caption: Binning
   :toctree:
   :nosignatures:

   ~westpa.BinMapper
   ~westpa.BinMapperBase
   ~westpa.RectilinearBinMapper
   ~westpa.VoronoiBinMapper
   ~westpa.MABBinMapper
   ~westpa.RecursiveBinMapper


Resampling
----------

.. autosummary::
   :caption: Resampling
   :toctree:
   :nosignatures:

   ~westpa.Resampler
   ~westpa.ResamplerBase
   ~westpa.HuberKimResampler
   ~westpa.MultinomialResampler
   ~westpa.ResidualResampler
   ~westpa.StratifiedResampler
   ~westpa.SystematicResampler

