To Do List
==========

- Use automatic archive of rebound:
  https://rebound.readthedocs.io/en/latest/ipython_examples/SimulationArchiveRestart/

- Create spherical shell type of object.

- Simulate different transit depths superposing spherical shells
  responding differently to wavelengths.

- Create an animation showing how the appearance of the object change
  as the observer, the star or the orientation of the body change [DONE].

- When purging points in Spangler calculate distances using the arc
  distance instead of the cartesian (euclidean) distance.

- Include Rossiter-McLaughlin effect, https://arxiv.org/pdf/1709.00680.pdf.

- Include radial velocities calculation.  For each object generate a
  curve of radial velocity in the direction of the observer.  Radial
  velocities could be calculated for individual spangles so we can
  simulate spectra modifications due for instance to stellar rotation.

- Make code Python PEP8 compliant:
  https://ellibrodepython.com/python-pep8#:~:text=La%20PEP8%20es%20una%20gu%C3%ADa,que%20una%20l%C3%ADnea%20debe%20tener.
  https://peps.python.org/pep-0008/

- Use standard docstring conventions:
  https://peps.python.org/pep-0257/
  https://numpydoc.readthedocs.io/en/latest/format.html
  https://www.datacamp.com/tutorial/docstrings-python

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

- For intersections, transit and occult use the shape of the figure
  for defining the region on which this intersections occur. Use
  ConvexHull [DONE].

- REBOUND integration [DONE].

- Animation command using animated gif instead of ffmpeg [DONE].

- Use snake_case naming convention [DONE]

- Allow to change dynamically the orientation of the ring. Create an
  interact notebook to show orientation of the ring [DONE].

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

- Use already store sampling sets for spangling objects [DONE]

Other Ideas
===========

- Table of spangles: after updating the positions of all spangles in
  the system, create a table (dataframe) with all the spangles, their
  positions and properties in order to use it for computing states.

- Lightning synthesis: after generating the table of spangles, update
  the lighting states of the spangles computing light blocking
  (shadow) from other spangles.

  - Starting from the body to which spangle belong, explore the
    spangles of the parent and the children of the body.  Not include
    other bodies in the system (they are too far). Assume a parallel
    shadow.

- Observer synthesis: after generating the table of spangles, update
  the states of the spangles with respect to an observer by checking
  which spangles are blocked by another spangles, ie. which are
  visible and which are invisible.

  - Compute rectangular coordinates of the spangles in the observer
    references frame.

  - Check if center of spangle is at a distance of any other spangle
    smaller than their radius.

- Atmosphere: create a spherical body around a given body whose
  spangles be semitransparent.  When the light passing through the
  semitransparent body, it absorbes.

