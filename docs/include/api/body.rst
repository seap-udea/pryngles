Body Module
========================

The ``pryngles.body`` module provides the base classes for create and assign physical, orbital and 
optical parameters to an astrophysical body within a simulation for an arbitrary extrasolar system (see :doc:`system`).

.. autoclass:: pryngles.body.Body
   :members:
   :member-order: bysource

Inheritance from Body Class
----------------------------------

Base class :data:`~ body.Body` provides fundamental attributes for all astrophysical bodies in a simulated extrasolar system. 
Subclasses :data:`~ body.Star`, :data:`~ body.Planet`, and :data:`~ body.Ring` inherit from :data:`~ body.Body` and specialize its functionality: 

.. autoclass:: pryngles.body.Star
   :member-order: bysource
   :members:

.. autoclass:: pryngles.body.Planet
   :member-order: bysource
   :members:

.. autoclass:: pryngles.body.Ring
   :member-order: bysource
   :members:

.. autoclass:: pryngles.body.Observer
   :member-order: bysource
   :members:

.. code-block:: python
   :caption: Examples

   >>> # Let's create the body components of an exoplanetary system
   >>> Star = pr.Star(kind = 'Star', defaults = pr.STAR_DEFAULTS, parent = None)
   >>> Planet = pr.Planet(kind = 'Planet', defaults = pr.PLANET_DEFAULTS, parent = Star)
   >>> Ring = pr.Ring(kind = 'Ring', defaults = pr.RING_DEFAULTS, parent = Planet)

