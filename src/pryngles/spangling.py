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

Spangling_doc="""A Spangling associated to an object or set of objects.
    
Initialization attributes:

    df : Pandas DataFrame: 
        Dataframe containing the information on the spangling.
        
    body_hash

Secondary attributes:

    hash: string
        Hash of the object, ie. a unique string identifying the object 
        (see hash Python function)

Public methods:

    update_body(props):
        Update a given property.
"""

#Columns of spangling
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
        "ns_equ":[0,0,1],
        "ns_ecl":[0,0,1],
        "ns_obs":[0,0,1],
        #Optical parameters
        "asp":0.0,
        "albedo_gray_normal":1,
        "tau_gray_optical":0,
        #Spangle state
        "unset":1, #State has not been set
        "visible":0, #The spangle is visible from observer
        "shadow":0, #The spangle is in the shadow of other spangle
        "illuminated":0, #The spangle is illuminated
        "transit":0, #The spangle is transiting
        "indirect":0, #The spangle is indirectly illuminated
        "occult":0, #The spangle is occulted by a light source
    }
)

class Spangling(object):
    
    def __init__(self,nspangles=0,body_hash=None):

        #Attributes
        self.nspangles=nspangles
        
        #Update default values
        self.defaults=deepcopy(SPANGLING_COLUMNS)
        
        if body_hash:
            self.defaults.update(dict(body_hash=body_hash))
        
        #Declare dataframe
        if self.nspangles>0:
            self.df=pd.DataFrame([list(self.defaults.values())]*nspangles,columns=self.defaults.keys())
        else:
            self.df=pd.DataFrame(columns=self.defaults.keys())

from IPython.display import display, HTML


# ### The end

Spangling.__doc__=Spangling_doc

