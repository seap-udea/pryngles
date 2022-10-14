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

# # Pryngles module: physics
# 
# Template of a module

from pryngles import *

# ## External modules

#@external
import numpy as np
import rebound as rb
from tqdm import tqdm
#@end:external

# ## Orbit of Hierarchical N-body system

# ### Docstring

Orbit_doc="""Hierarchical N-body system.

Initialization parameters+:

    m1,m2: mixed, default = 1:
        If m1 (m2) is a number, this is the mass of the first body.
        If m1 (m2) is an instance of Orbit, this is a children system.

    elements: dictionary:
        Dictionary with the elements provided to indicate the location of the relative 
        vector of the system.
        
        Valid orbital elements are:
            
            a=1,f=0,e=0: 
                Semi major axis, true anomaly, eccentricity.
            
            omega=0,inc=0,Omega=0:
                Periapsis argument, inclination, longitude of the ascending node.

            M=0,E=0,T=0:
                Mean anomaly, eccentric anomaly, time of periapsis passage.
                
            theta=0,l=0:
                True longitude (Omega + omega + f), mean longitude (Omega + omega + M).

    R: array (3), default = [0,0,0]
        Initial position of the center of mass.
    
    V: array (3), default = [0,0,0]
        Initial velocity of the center of mass.
    
Secondary attributes:

    m1, m2: float:
        Masses of the particles.

    Mtot: float:
        Total mass of the system.
    
    sub_sim: Rebound Simulation:
        Simulation corresponding to the motion of the two particles.
        
    sim: Rebound simulation:
        Simulation corresponding to all particles in the system.

Key methods:

    get_positions():
    
        Returns:
            states: list (N) of dictionaries:
                List of dictionaries having the state vector x,y,z,vx,vy,vz of 
                each particle.

All methods:
    See Misc.get_methods(Orbit)
    
Examples:

    Complex system: two binaries orbited by a further away object:
    
        S1=Orbit(m1=1,m2=1,a=1,e=0.7,M=0)
        S2=Orbit(m1=1,m2=1,a=1,e=0,M=0)
        S3=Orbit(S1,S2,a=5,e=0)
        S4=Orbit(S3,m2=1,a=20,e=0,E=45*Consts.deg)
        S4.ensamble_system()
        Plot.animate_rebound(S4.sim)
        
    Planet with a moon:
        units=["au","msun","yr"]
        orb=Orbit(
            m1=1,
            m2=Orbit(m1=1e-3,m2=1e-7,a=0.5,e=0.0,units=units),
            units=units,
            a=20,e=0.0)
        orb.calculate_orbit()
        sim,states=orb.get_states()
        Plot.animate_rebound(sim)
        
    Simple system:
        units=["au","msun","yr"]
        sim,states=Orbit(m1=1.0,m2=1e-3,units=units,a=0.1,e=0.0).calculate_orbit().get_states()
        
"""

# ### Class HNbody

class Orbit(PrynglesCommon):
    
    def __init__(self,m1=1,m2=1,
                 R=np.array([0,0,0]),V=np.array([0,0,0]),
                 units=None,
                 **elements):
        
        #Global simularion
        self.sim=rb.Simulation()
        if units:
            self.units=units
        else:
            self.units=["au","msun","yr2pi"]
        self.sim.units=self.units

        #Particles
        self.p1=m1
        self.p2=m2
        
        #Periods
        self.Ps=[]
        
        #Check first system
        qmixed=False
        if isinstance(self.p1,Orbit):
            self.m1=self.p1.Mtot
            qmixed=True
        elif isinstance(self.p1,float) or isinstance(self.p1,int):
            self.m1=self.p1
        else:
            raise ValueError(f"Type of first componente ({type(m1)}) not recognized.  It should be a float or an Orbit instance.")
        
        #Check second system
        if isinstance(self.p2,Orbit):
            self.m2=self.p2.Mtot
            qmixed=True
        elif isinstance(self.p2,float) or isinstance(self.p2,int):
            self.m2=self.p2
        else:
            raise ValueError(f"Type of Second component ({type(m2)}) not recognized.  It should be a float or an Orbit instance.")
                
        if not qmixed and (sum(R)!=0 or sum(V)!=0):
            raise ValueError(f"You cannot provide a center of mass position and velocity for a non-mixed system.")
        
        #Total mass
        self.Mtot=self.m1+self.m2
        
        #Add initial elements to attributes
        for element in elements:
            if element in REBOUND_ORBITAL_PROPERTIES:
                self.__dict__[element]=elements[element]
            else:
                raise ValueError(f"Element {element} not identified.")
                
        #Update states
        self._update_states(R,V)
        
    def _update_states(self,R=np.array([0,0,0]),V=np.array([0,0,0])):
        """Update state of the particles
        """        
        #Create rebound obtions
        self._rb_options={k:v for k,v in self.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}  
        self._rb_options.update(dict(m=0))
        
        #Create rebound simulation
        sim=rb.Simulation()
        sim.units=self.units
        sim.add(m=self.Mtot)
        sim.add(**self._rb_options)
        
        #Relative vector
        self.r=np.array(sim.particles[1].xyz)
        self.v=np.array(sim.particles[1].vxyz)
        
        del sim
        
        #Calculate positions of components
        self.r1=R-self.m2/self.Mtot*self.r
        self.v1=V-self.m2/self.Mtot*self.v
        
        if isinstance(self.p1,Orbit):
            self.p1._update_states(self.r1,self.v1)
            
        self.r2=R+self.m1/self.Mtot*self.r
        self.v2=V+self.m1/self.Mtot*self.v                
        if isinstance(self.p2,Orbit):
            self.p2._update_states(self.r2,self.v2)
            
        #Create a simulation of this system
        self.sub_sim=rb.Simulation()
        self.sub_sim.units=self.units
        self.sub_sim.add(m=self.m1,
                         x=self.r1[0],y=self.r1[1],z=self.r1[2],
                         vx=self.v1[0],vy=self.v1[1],vz=self.v1[2])
        self.sub_sim.add(m=self.m2,
                         x=self.r2[0],y=self.r2[1],z=self.r2[2],
                         vx=self.v2[0],vy=self.v2[1],vz=self.v2[2])

    def calculate_orbit(self,sim=None):
        """Ensamble Hierarchical N-body system.
        
        Parameters:
            sim: Rebound Simulation, default = None:
                This is used for recursion purposes.  
                Normally is set to None.
        
        Return:
            orbit: Orbit object:
                Used for nested purposes, 
                Example: orbit.calculate_orbit().get_states()
        """
        if sim is None:
            sim=self.sim
            
        if isinstance(self.p1,Orbit):
            self.p1.calculate_orbit(sim)
        else:
            sim.add(m=self.m1,
                    x=self.r1[0],y=self.r1[1],z=self.r1[2],
                    vx=self.v1[0],vy=self.v1[1],vz=self.v1[2])
            
        if isinstance(self.p2,Orbit):
            self.p2.calculate_orbit(sim)
        else:
            sim.add(m=self.m2,
                    x=self.r2[0],y=self.r2[1],z=self.r2[2],
                    vx=self.v2[0],vy=self.v2[1],vz=self.v2[2])
            
        for p in sim.particles[1:]:
            self.Ps+=[p.P]
        return self
            
    def get_states(self):
        """Get positions of particles in the system
        
        Returns:
            states: list (N) of dictionaries:
                List of dictionaries having the state vector x,y,z,vx,vy,vz of 
                each particle.        
        """
        states=[]
        for p in self.sim.particles:
            states+=[
                dict(m=p.m,x=p.x,y=p.y,z=p.z,vx=p.vx,vy=p.vy,vz=p.vz)
            ]
        return self.sim,states

#@test:template
#@end

