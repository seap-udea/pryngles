Orbit Module
========================

This is our Orbital Mechanics Module. This module provides classes and utilities for building and simulating hierarchical N-body systems, as this is the most common framework to model extrasolar systems. 
It includes functionality for defining orbital bodies, constructing orbits, and performing celestial mechanics calculations. 

.. note::
   ``Pryngles`` implements the facilities of the well known `REBOUND <https://rebound.readthedocs.io/en/latest/>`_ package to optimize its dynamical solutions. 

.. important::
   By hierarchical N-body systems, we mean systems whose dynamics can be represented by the 2-body subsystem clustering approximation.  
   The next figure may illustrate you the kind of systems we want to reproduce.

   .. image:: images/hierarchical.png
      :align: center
      :scale: 80 %

   Make a smart construction of your system, i.e., build your system by segmenting it into nested 2-body subsystems.  
   Please follow the next scheme convention for organize your hierarchical system

   .. image:: images/TYC_diagram.png
      :scale: 50 %
      :align: center 

.. autoclass:: pryngles.orbit.OrbitUtil
   :members:
   :member-order: bysource


.. autoclass:: pryngles.orbit.Orbit
   :members:
   :member-order: bysource

