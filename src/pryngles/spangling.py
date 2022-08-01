##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: < template >

from pryngles import *

# ## External modules

# ## The body class
# 
# The Body class is one of the most important classes in the package. 

Class_doc="""A general body.  This calss is not intended to be used independently, just for inheritance purposes.
    
Initialization attributes:

    kind : string
        One of the kind of bodies defined in the package (see _BODY_KINDS)

Secondary attributes:

    hash: string
        Hash of the object, ie. a unique string identifying the object 
        (see hash Python function)

Public methods:

    update_body(props):
        Update a given property.
"""

class Class(object):
    
    def __init__(self):
        pass

Class.__doc__=Class_doc

# ## Testing


