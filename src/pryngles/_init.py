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

# # Pryngles module

# ## Module of basic functions
# 
# Goals of the module:
# - Future `__init__`

##HEADER
from pryngles import *

import sys
IN_JUPYTER='ipykernel' in sys.modules

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

# ## Package macros use across the modules

from collections import OrderedDict as odict
import os
import sys
import numpy as np
import pkg_resources
from rebound import units

# ## Global constants
# 
# All global constants in Pryngles have a capital name and a lowercase version in the class util.

#Root directory
try:
    __rootdir__=os.path.abspath(os.path.dirname(__file__))
except:
    import IPython
    __rootdir__=os.path.abspath('')

#Version
__version__=pkg_resources.require('pryngles')[0].version

# ### Mathematical and physical constants

class Consts(object):pass

#Constants and units
#Mathematical constants
RAD=180/np.pi
DEG=1/RAD
Consts.rad=180/np.pi
Consts.deg=1/Consts.rad
Consts.ppm=1e6 #parts per million
Consts.ppb=1e9 #parts per billion

#Physical constants
GSI=units.convert_G(["m","s","kg"]) # G constant in SI units

from rebound import units
for const in "times","lengths","masses":
    values=eval(f"units.{const}_SI.copy()")
    for key in values:
        exec(f"Consts.{key}=values[key]")

#Size of reference objects 
Consts.rearth=6378.137e3 #m, volumetric mean radius
Consts.rsun=695700e3 #m, nominal solar radius
Consts.rjupiter=71492e3 #m, equatorial radius
Consts.rsaturn=60268e3 #m

# ## Clases

util_doc="""
Util routines.

This is a set of util routines intended for a diversity of purposes.

Routines included:

    get_data(file)
""";

class util(object):
    def get_data(path):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            datafile: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(__rootdir__,'data',path);
util.__doc__=util_doc

class PrynglesCommon(object):
    def __str__(self):
        return str(self.__dict__)    

