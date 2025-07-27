.. _system:

System Module
=============

This module provides the `System` class for creating and managing planetary systems within the **Pryngles** package. The `System` class integrates with the `rebound` library to simulate N-body dynamics and uses :any:`spangler.Spangler` discretizing method to model the surfaces of astrophysical bodies, enabling calculations of light scattering, transit, shadowing, and visibility.

.. note::
   In ``Pryngles``, the `System` class orchestrates the simulation of multiple bodies (:any:`body.Star`, :any:`body.Planet`, :any:`body.Ring`) and their interactions, combining orbital dynamics with radiative transfer calculations for astrophysical modeling.

.. autoclass:: pryngles.system.System
   :members:
   :member-order: bysource