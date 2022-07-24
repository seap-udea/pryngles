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

# ## External modules

# ## Ring Class

class RingDefaults(object):
    """
    These are the default attributes for bodies of the kind 'Ring'.
    
    DEVELOPER:
        You may add new attributes as the model gets more complex.
        Please document properly each attribute.
        
    orbit:
    
        (In current version, Ring body does not have orbit attributes)
        
    physics:
    
        fi: float [adimensional], default = 1
            Fraction of the radius of the primary object where ring stars.
            
        fe: float [adimensional], default = 1
            Fraction of the radius of the primary object where ring ends.
            
        i: float [radians], default = 0
            Inclination of the ring with respect to the equator of the primary
            object.
            
    optics:
    
        nspangles: int, default = 1000
            Number of spangles on which the object will be discretized.
    """
    orbit=dict()
    
    physics=dict(fi=1.0,fe=2.0,i=0.0)
    
    optics=dict(nspangles=1000)

BODY_KINDS+=["Ring"]

class Ring(Body):
    """Class Planet.
    
    See Body class documentation.
    
    Additional public attributes:
    
        physics.ri: float [rebound length unit]
            Radius of the inner border of the ring

        physics.re: float [rebound length unit]
            Radius of the outer border of the ring
    
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
        self.primary=primary
        self.primary._update_childs(self)
        self._update_parent(self.primary)
        
        #Update properties
        self.update_body(**self.__dict__)
        
    def update_body(self,**pars):
        Body.update_body(self,**pars)
        
        #Here place the commands to update this kind of body
        self.physics.ri=self.physics.fi*self.primary.physics.radius
        self.physics.re=self.physics.fe*self.primary.physics.radius


