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

# # Pryngles module: Constants and Enumerators
# 
# This module contains all the constants required by other modules in the package.

from pryngles import *

# ## External modules

import numpy as np
from rebound import units

# ## Constants and enumerators by Module

# ### System related

import os
#Root directory
try:
    FILE=__file__
    ROOTDIR=os.path.abspath(os.path.dirname(FILE))
except:
    import IPython
    FILE=""
    ROOTDIR=os.path.abspath('')
    
IN_JUPYTER='ipykernel' in sys.modules

# ### Body related

BODY_KINDS=[]

# ### Spangle related

#Type of spangle
SOLID_SPANGLE=0
GRANULAR_SPANGLE=1
ATMOSPHERIC_SPANGLE=2

#Reference system
EQU=0
ECL=1

# ## Constants

class Consts(object):pass

#Mathematical constants
Consts.rad=180/np.pi
Consts.deg=1/Consts.rad
Consts.ppm=1e6 #parts per million
Consts.ppb=1e9 #parts per billion

#Physical constants
GSI=units.convert_G(["m","s","kg"]) # G constant in SI units
for const in "times","lengths","masses":
    values=eval(f"units.{const}_SI.copy()")
    for key in values:
        exec(f"Consts.{key}=values[key]")

#Size of reference objects
Consts.rearth=6378.137e3 #m, volumetric mean radius
Consts.rsun=695700e3 #m, nominal solar radius
Consts.rjupiter=71492e3 #m, equatorial radius
Consts.rsaturn=60268e3 #m


