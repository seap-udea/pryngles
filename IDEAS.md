Ideas
=====

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

