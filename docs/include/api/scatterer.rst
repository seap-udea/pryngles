Scatterer Module
========================

This module is intended for implementing light scattering models in exoplanetary atmospheres and surfaces, with the goal of reproducing their photometric signatures.

.. note::

   As you may know, ``pryngles`` works based on the idea of sampling a surface with several discrete units that we named `Spangles`. For each `Spangle`, we define its `Scatterer` object for photometry modelling purposes


.. seealso::
   
   For the science behind the definition and computation of our methods, please refer to the following literature

   | **[1]** Russell, H. N. (1916). On the Albedo of Planets and Their Satellites. `10.1086/142244 <https://ui.adsabs.harvard.edu/abs/1916ApJ....43..173R/abstract>`_
   | **[2]** Sobolev, V. V. (1975). Light Scattering in Planetary Atmospheres.  
   | **[3]** Yanovitskii, E. G. (1973). Spherical Albedo of a Planetary Atmosphere. `1973SvA....16..687Y <https://ui.adsabs.harvard.edu/abs/1973SvA....16..687Y/abstract>`_

   For an explanation of our implementation, please refer to our following paper

   **[1]** Zuluaga, J. I., Sucerquia, M., & Alvarado-Montes, J. A. (2022). 
   `The bright side of the light curve: A general photometric model of non-transiting exorings`. 
   Astronomy and Computing 40 (2022) 100623. `arXiv:2207.08636 <https://arxiv.org/abs/2207.08636>`_

.. autoclass:: pryngles.scatterer.Scatterer
   :exclude-members:

.. autoclass:: pryngles.scatterer.NeutralSurface
   :members:

.. autoclass:: pryngles.scatterer.BlackBodySurface
   :members:

.. autoclass:: pryngles.scatterer.LambertianGraySurface
   :members:

.. autoclass:: pryngles.scatterer.LambertianGrayAtmosphere
   :members:



.. .. automodule:: pryngles.scatterer
..    :members:
..    :undoc-members:
..    :show-inheritance: