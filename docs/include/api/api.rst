API Reference
========================

| ``Pryngles`` is a highly modular package divided into modules and data files.  
| In designing the package, we have followed the convention of creating, for each module, a corresponding ``Class`` with the same name as the module. 
Thus, for instance, the module :doc:`Spangler <spangler>` contains a general-purpose class :any:`spangler.Spangler` that includes routines intended for sampling an astrophysical body surface.

| ``Pryngles`` has two types of modules:

- **General modules**: These modules contain *static classes*, i.e. classes that cannot be instantiated and their methods work as regular functions.

  - :doc:`consts`: Constants and other methods associated with them.
  - :doc:`Science <science>`/:doc:`scatterer`: Classes and methods involving useful scientific calculations.

- **Core modules**: These are the key modules in the package. They contain regular classes and perform the most specialized functions in the package.

  - :doc:`sampler`: Module containing the classes and methods to generate samples of points.
  - :doc:`spangler`: Classes and methods to spangle bodies.
  - :doc:`body`: Classes and methods used to create and spangle astrophysical bodies.
  - :doc:`orbit`: Main module to create hierarchical N-body system.
  - :doc:`system`: This is the most important module. It contains the classes and methods used to define a planetary system, spangle it, and calculate light-curves.
  .. - :doc:`legacy`: This module contains all the classes and routines of the legacy ``RingedPlanet`` interface to Pryngles. This interface will be deprecated in version 2.0.

.. toctree::
   :titlesonly:
   :hidden:

   body
   consts
   orbit
   sampler
   scatterer
   science
   spangler
   system
