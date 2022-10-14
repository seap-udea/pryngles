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
import re

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

# ### Planet related

#Default colors
PLANET_COLOR="#D3CAAF"
STAR_COLOR="#FDFE7D"
RING_COLOR="#DBD5D1"

# ## System related

REBOUND_ORBITAL_PROPERTIES=dict(
    #Mass
    m=0,
    #Cartesian coordinates
    x=0,y=0,z=0,vx=0,vy=0,vz=0,
    #Semi major axis, true anomaly, eccentricity
    a=1,f=0,e=0,
    #Periapsis argument, inclination, longitude of the ascending node
    omega=0,inc=0,Omega=0,
    #Mean anomaly, eccentric anomaly, time of periapsis passage
    M=0,E=0,T=0,
    #true longitude (Omega + omega + f), mean anomaly (Omega + omega + M)
    theta=0,l=0,
)

# ## Constants class
# 
# This class contains several physical and astronomical constants.  It takes constants from other packages (rebound, astropy) and define other constants based on reliable sources.

class Consts(object):pass

#Mathematical constants
Consts.rad=180/np.pi
Consts.deg=1/Consts.rad
Consts.ppm=1e6 #parts per million factor
Consts.ppb=1e9 #parts per billion factor

#Physical constants
GSI=units.convert_G(["m","s","kg"]) # G constant in SI units
for const in "times","lengths","masses":
    values=eval(f"units.{const}_SI.copy()")
    for key in values:
        exec(f"Consts.{key}=values[key]")

#Size of reference objects
Consts.rearth=6378.137e3 #m, volumetric mean radius, source: 
Consts.rsun=695700e3 #m, nominal solar radius, source: 
Consts.rjupiter=71492e3 #m, equatorial radius, source: 
Consts.rsaturn=60268e3 #m, equatorial radius, source: 

#For compatibility purposes with legacy: remove when legacy is retired
RAD=Consts.rad
DEG=Consts.deg

def get_physical():
    import pryngles as pr
    all_constants=[]
    for key in Consts.__dict__.keys():
        patterns = "^[a-z]+$"
        if re.search(patterns,key):
            all_constants+=[key]
    return sorted(all_constants)

def get_all():
    import pryngles as pr
    all_constants=[]
    for key in pr.__dict__.keys():
        patterns = "^[A-Z_]+$"
        if re.search(patterns,key):
            all_constants+=[key]
    return sorted(all_constants)

Consts.get_physical=get_physical
Consts.get_all=get_all


