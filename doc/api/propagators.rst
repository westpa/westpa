Propagators
===========

.. autoclass:: westpa.Propagator()
   :members: __call__

.. autoclass:: westpa.SerialPropagator
   :members: propagate, make_segment_dir

.. autoclass:: westpa.VectorizedPropagator
   :members: propagate, make_segment_dir

.. autoclass:: westpa.AmberPropagator

.. autoclass:: westpa.GROMACSPropagator

.. autoclass:: westpa.OpenMMPropagator
   :members: add_reporter
