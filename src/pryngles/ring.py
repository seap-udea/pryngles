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

# # Pryngles module: rings

from pryngles import *
sci=Science

# ## External modules

import spiceypy as spy
import math as mh
import numpy as np
import copy

# ## Ring default properties

class RingDefaults(object):
    """
    These are the default attributes for bodies of the kind 'Ring'.
    
    DEVELOPER:
        You may add new attributes as the model gets more complex.
        Please document properly each attribute.
        
    orbit:
    
        (In current version, Ring body does not have orbit attributes)
        
    physics:
    
        fi: float [adimensional], default = 1:
            Fraction of the radius of the primary object where ring stars.
            
        fe: float [adimensional], default = 1:
            Fraction of the radius of the primary object where ring ends.
            
        i: float [radians], default = 0:
            Inclination of the ring with respect to the equator of the primary
            object.
            
        roll: float [radians], default = 0:
            Roll angle.  This is the angle with respect to ecliptic x-axis in which 
            the normal to the ring plane is rotated.
            
    optics:
    
        nspangles: int, default = 1000: 
            Number of spangles on which the object will be discretized.
            
        albedo_gray_normal: float. default = 1: 
            Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.
            
        tau_gray_optical: float. default = 0:
            Gray (wavelength indpendent) Optical depth of the spangle.  
            If 0 the spangle is entirely opaque to all wavelength, despite its type.            
    """
    orbit=dict()
    
    physics=dict(fi=1.0,fe=2.0,i=0.0,roll=0.0)
    
    optics=dict(nspangles=1000,
                albedo_gray_normal=1,
                tau_gray_optical=0
               )

# ## Ring Class

BODY_KINDS+=["Ring"]

class Ring(Body):
    """Class Planet.
    
    See Body class documentation.
    
    Additional public attributes:
    
        ri: float [rlu]:
            Radius of the inner border of the ring

        re: float [rlu]:
            Radius of the outer border of the ring
            
        spangles: list of Spangle objects:
            List of spangles covering the surface (for spangle options see Spangle documentation)
    
    Override methods:
    
        update_body(**pars):
            This method compute additional attributes like (see above).
    """
    def __init__(self,
                 primary=None,
                 orbit=RingDefaults.orbit,
                 physics=RingDefaults.physics,
                 optics=RingDefaults.optics
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,RingDefaults,"Ring",primary,orbit,physics,optics)
        
        #Check primary
        if self.primary is None:
            raise ValueError(f"Primary not provided and mandatory for {self.kind}.")
        #self.primary=primary
        #self.primary._update_childs(self)
        #self._update_parent(self.primary)
        
        #Update properties
        self.update_body(**self.__dict__)
        
    def _set_transform(self):
        """
        Set transformation matrices
        """
        #Normal vector to ring
        self.nr_ecl=np.array([0,mh.sin(self.physics.i),mh.cos(self.physics.i)]) 
        self.nr_equ=np.array([0,0,1]) 
        
        #Build transformation matrices
        ex_equ=np.array([1,0,0]) 
        ey_equ=np.array([0,mh.cos(self.physics.i),-mh.sin(self.physics.i)])
        ez_equ=self.nr_ecl
        self.M_equ2ecl=np.array(list(np.vstack((ex_equ,ey_equ,ez_equ)).transpose())).reshape((3,3))
        self.M_ecl2equ=spy.invert(self.M_equ2ecl)


# ### Update body

def update_body(self,**pars):
    Body.update_body(self,**pars)

    #Here place the commands to update this kind of body
    self.ri=self.physics.fi*self.primary.physics.radius
    self.re=self.physics.fe*self.primary.physics.radius

    #Set transformation matrices
    self._set_transform()

Ring.update_body=update_body


# ### Spangle body
# 
# This method creates the spangles taht will cover the body.

def spangle_body(self,seed=0):
    #Create spangler
    self.sp=Spangler(N=self.optics.nspangles)
    
    #Limits of the ring (normalized to re)
    uri=self.ri/self.re
    ure=1
    
    #Generate spangles
    self.sp.gen_ring([
        [0.0,uri],#Internal gap
    ],boundary=0)
    self.sp.purge_sample()
    
    #Body spangles
    self.spangles=np.array([],dtype=Spangle)
    for i in range(self.sp.N):

        #Create spangle
        spangle=Spangle()
        
        #Type of spangles
        spangle.type=GRANULAR_SPANGLE
        
        #Coordinates of the spangle
        xyz_equ=np.append(self.sp.ss[i]*self.re,0) #Complete coordinates with z = 0
        xyz_ecl=spy.mxv(self.M_equ2ecl,xyz_equ)
        
        rqf_equ=np.append(self.sp.pp[i],0) #Complete coordinates with phi = 0
        rqf_equ[0]*=self.re
        rqf_ecl=sci.xyz2rqf(xyz_ecl) #Convert equatorial to spherical
        spangle.set_position([xyz_equ,xyz_ecl],[rqf_equ,rqf_ecl])

        #Unitary vector normal to the spangle
        ns_equ=np.array([0,0,1])
        ns_ecl=spy.mxv(self.M_equ2ecl,ns_equ)
        spangle.set_orientation([ns_equ,ns_ecl])
            
        #Spangle Area
        spangle.asp=self.sp.aes*self.re**2
        
        #Optical properties
        spangle.set_optical(
            albedo_gray_normal=self.optics.albedo_gray_normal,
            tau_gray_optical=self.optics.tau_gray_optical
        )
        
        self.spangles=np.append(self.spangles,copy.deepcopy(spangle))
        del spangle

Ring.spangle_body=spangle_body


