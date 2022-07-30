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

# # Pryngles module: Spangle

from pryngles import *

# ## External modules

import numpy as np

# ## Enumerators and constantes

# ## Spangle class

Spangle_doc=f"""    Spangles are the basic unit of Pryngles.  This class is initially intended to contain 
    most of the relevant physics of light scattering and polarization.
    
    Initialization attributes:
        
    Basic attributes:
    
        type: integer. default = SOLID_SPANGLE :
            Type of spangle (SOLID_SPANGLE, GRANULAR_SPANGLE, ATMOSPHERIC_SPANGLE, see consts module).
    
        xyz: list or array (2).  default = [[0,0,0]]*2: 
            Cartesian coordinates of the spangle in the simulation.  
            The 0-component (EQU component) are the coordinates in the equatiorial system
            of the object to which the spangle belong.
            The 1-component (ECL component) are the coordinates in the ecliptic system (orbit)
            of the object to which the spangle belong.
            This class is not responsible for changing from one system of coordinates to the 
            other.
            
        rqf: list or array (2).  default = [[0,0,0]]*2: 
            Spherical coordinates of the spangle in the simulation (see xyz for details).  
        
        ns: list of arrays (2). default = [[0,0,1]]*2:
            Unitary vector normal to the surface of the spangle in the EQU and ECL system (see xyz for details).
        
    Physical attributes:
    
        albedo_gray_normal: float. default = 1: 
            Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.
            
        tau_gray_optical: float. default = 0:
            Gray (wavelength indpendent) Optical depth of the spangle.  
            If 0 the spangle is entirely opaque to all wavelength, despite its type.
"""

class Spangle(PrynglesCommon):
    def __init__(self):

        #Basic
        self.type=SOLID_SPANGLE
        
        #Position 
        self.xyz=[np.array([0,0,0])]*2 #Cartesian coordinates
        self.rqf=[np.array([0,0,0])]*2 #Spherical or cylindrical coordinates
        
        #Area
        self.asp=1
        
        #Orientation
        ns=[np.array([0,0,1])]*2
        
        #Completely reflective spangle
        self.albedo_gray_normal=1.0
        
        #Completely opaque spangle (despite type)
        self.tau_gray_optical=0.0
        
    def set_position(self,xyz,rqf):
        
        try:
            len(xyz)
        except:
            raise IndexError(f"The cartesian coordinates of the spangle should be a list.  You provided {xyz}")

        if len(xyz)!=2:
            raise IndexError(f"The list of cartesian coordinates should be 2.  You provided {xyz}")
            
        if len(xyz[EQU])!=3 or len(xyz[ECL])!=3:
            raise IndexError(f"The cartesian coordinates of the spangle should have len 3.  You provided {xyz}")

        try:
            len(rqf)
        except:
            raise IndexError(f"The spherical coordinates of the spangle should be a list.  You provided {rqf}")
        
        if len(rqf)!=2:
            raise IndexError(f"The list of spherical coordinates should be 2.  You provided {rqf}")

        if len(rqf[EQU])!=3 or len(rqf[ECL])!=3:
            raise IndexError(f"The spherical coordinates of the spangle should have len 3.  You provided {rqf}")

        self.xyz=xyz
        self.rqf=xyz
    
    def set_orientation(self,ns):
        try:
            len(ns)
        except:
            raise IndexError(f"The normal vector to the spangle should be a list.  You provided {ns}")

        if len(ns[EQU])!=3 or len(ns[ECL])!=3:
            raise IndexError(f"The normal vector to the spangle should have len 3.  You provided {ns}")
                
        self.ns=ns
        
    def set_optical(self,**props):

        for prop in props.keys():
            if prop not in self.__dict__:
                raise KeyError(f"Optical property '{prop}' not recognized in this class")

        prop="albedo_gray_normal"
        if prop in props:
            if props[prop]>1 or props[prop]<0:
                raise ValueError(f"{prop} cannot be larger than 1 or smaller than 0.  You provided {props[prop]}")
            self.__dict__[prop]=props[prop]
        
        prop="tau_gray_optical"
        if prop in props:
            if props[prop]<0:
                raise ValueError(f"{prop} cannot be smaller than 0.  You provided {props[prop]}")
            self.__dict__[prop]=props[prop]

Spangle.__doc__=Spangle_doc


