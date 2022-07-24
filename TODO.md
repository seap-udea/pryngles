To Do List
==========

- Create a full test notebook.

- Create a version file.

- Include the ACL code: https://ascl.net/code/v/3289 (when accepted).

- Animation command using animated gif instead of ffmpeg.

- Make code Python PEP8 compliant:
  https://ellibrodepython.com/python-pep8#:~:text=La%20PEP8%20es%20una%20gu%C3%ADa,que%20una%20l%C3%ADnea%20debe%20tener.
  https://peps.python.org/pep-0008/

- Use snake_case naming convention.

- Use standard docstring conventions:
  https://peps.python.org/pep-0257/
  https://numpydoc.readthedocs.io/en/latest/format.html
  https://www.datacamp.com/tutorial/docstrings-python

- Allow to change dynamically the orientation of the ring.

- Create an interact notebook to show orientation of the ring.

- Include units transformation inside system class for providing
  properties.

- Include photomeric filters, http://svo2.cab.inta-csic.es/theory/fps/

- Generate actual light curves including detector properties.

- Create a spangle class to describe the properties of each spangle.

- Create routines for creating the whole light-curve.

- Create interfaces for the components of the package that can be
  modified by the user: physics of scattering, physics of
  transmission.

  Example. scattering(lambda,eta,Ap): indicates the multiplying factors
  that apply on the amount of light received from direction lambda,
  outgoing in direction eta and at azimuth Ap.

- Incorporate wavelength dependent limb-darkening coefficients.

- Include forward scattering.

- Plot returns information to render a photorealistic image of the
  planet.  Returns the value of each point in a plot.

- Include the possibility of non-ringed planets (Nr=0).

- Include a citation command.

Future version [1.0]
====================

- [PARTIALLY DONE] REBOUND integration.

Done
====

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

