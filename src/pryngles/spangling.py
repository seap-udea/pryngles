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

# # Pryngles module: Spangling

# This module contains all the physics of light scattered on spangles

from pryngles import *

# ## External modules

import pandas as pd
from collections import OrderedDict as odict
from copy import deepcopy

# ## The Spangling class
# 
# This class contains a family of routines useful for spangling

Spangling_doc="""A general body.  This calss is not intended to be used independently, just for inheritance purposes.
    
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

SPANGLING_COLUMNS=odict(
    {
        "body_hash":"",
        #Type of spangle
        "type":SOLID_SPANGLE,
        #Cartesian coordinates
        "x_equ":0,"y_equ":0,"z_equ":0,
        "x_ecl":0,"y_ecl":0,"z_ecl":0,
        "x_obs":0,"y_obs":0,"z_obs":0,
        #Spherical coordinates
        "r_equ":0,"t_equ":0,"f_equ":0,
        "r_ecl":0,"t_ecl":0,"f_ecl":0,
        "r_obs":0,"t_obs":0,"f_obs":0,
        #Normal to spangle
        "nr_equ":[0,0,1],
        "nr_ecl":[0,0,1],
        "nr_obs":[0,0,1],
        #Optical constants
        "albedo_gray_normal":1,
        "tau_gray_optical":0
    }
)

class Spangling(object):
    
    def __init__(self,nspangles=0,body_hash=None):

        #Attributes
        self.nspangles=nspangles
        self.body_hash=body_hash
        
        #Update default values
        self.defaults=deepcopy(SPANGLING_COLUMNS)
        
        if body_hash:
            self.defaults.update(dict(body_hash=body_hash))
        
        #Declare dataframe
        if nspangles>0:
            self.df=pd.DataFrame([list(self.defaults.values())]*nspangles,columns=self.defaults.keys())
        else:
            self.df=pd.DataFrame(columns=self.defaults.keys())

from IPython.display import display, HTML


# ### The end

Spangling.__doc__=Spangling_doc

