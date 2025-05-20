Sampler Module
========================

This is our module to produce Fibonacci sampling over disk and sphere geometrical surfaces. This is the core of our :doc:`spangler` module, 
used to discretize the surface of astrophysical bodies and model their photometry

.. note::
   This module and the following class are based on the ``fibpy`` library by Martin Roberts.
   **Source code:** https://github.com/matt77hias/fibpy

.. autoclass:: pryngles.sampler.Sampler
   :members: plot, gen_circle, gen_ring, gen_sphere
   :member-order: bysource