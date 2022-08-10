To Do List
==========

- Use already store sampling sets for spangling objects.

- For intersections, transit and occult use the shape of the figure
  for defining the region on which this intersections occur. Use
  ConvexHull.

- Create spherical shell type of object.

- Simulate different transit depths superposing spherical shells
  responding differently to wavelengths.

- REBOUND integration [PARTIALLY DONE].

- Create an animation showing how the appearance of the object change
  as the observer, the star or the orientation of the body change.

- When purging points in Spangler calculate distances using the arc
  distance instead of the cartesian (euclidean) distance.

- Include Rossiter-McLaughlin effect, https://arxiv.org/pdf/1709.00680.pdf.

- Include radial velocities calculation.  For each object generate a
  curve of radial velocity in the direction of the observer.  Radial
  velocities could be calculated for individual spangles so we can
  simulate spectra modifications due for instance to stellar rotation.

- Animation command using animated gif instead of ffmpeg.

- Make code Python PEP8 compliant:
  https://ellibrodepython.com/python-pep8#:~:text=La%20PEP8%20es%20una%20gu%C3%ADa,que%20una%20l%C3%ADnea%20debe%20tener.
  https://peps.python.org/pep-0008/

- Use snake_case naming convention.

- Use standard docstring conventions:
  https://peps.python.org/pep-0257/
  https://numpydoc.readthedocs.io/en/latest/format.html
  https://www.datacamp.com/tutorial/docstrings-python

- Allow to change dynamically the orientation of the ring. Create an
  interact notebook to show orientation of the ring.

- Include photomeric filters, http://svo2.cab.inta-csic.es/theory/fps/

- Generate actual light curves including detector properties.

- Create routines for creating the whole light-curve.

- Create interfaces for the components of the package that can be
  modified by the user: physics of scattering, physics of
  transmission.

  Example. scattering(lambda,eta,Ap): indicates the multiplying factors
  that apply on the amount of light received from direction lambda,
  outgoing in direction eta and at azimuth Ap.

- Incorporate wavelength dependent limb-darkening coefficients.

- Include forward scattering for rings.

- Plot returns information to render a photorealistic image of the
  planet.  Returns the value of each point in a plot.

- Include the possibility of non-ringed planets (Nr=0).

- Include a citation command.

Done
====

- Create a spangle class to describe the properties of each spangle
  [DONE].

- Include units transformation inside system class for providing
  properties [DONE].

- Create a full test notebook [DONE].

- Create a version file [DONE]. 

- Include the ACL code: https://ascl.net/code/v/3289 (when accepted) [DONE].

- [DONE] Remove: Javascript Error: IPython is not defined.

- [DONE] Add __version__ variable to package.

- [DONE] Generate a water mark for the plots.
  Sol. Extra.prynglesMark()

- [DONE] Draw logo.
  Sol. Extra.drawPryngles()

- [DONE] Create the PyPI package.  See:
  https://packaging.python.org/en/latest/tutorials/packaging-projects/

  For the developers:

  - To build a distribution:

    python3 -m build
  
  - To upload the distribution:
  
    python -m twine upload --repository testpypi dist/* --verbose

- [DONE] Disable showing figure when plotting with showfig=0 in Google Colab.
