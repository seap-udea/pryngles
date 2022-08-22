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

# # Pryngles module: System

from pryngles import *

# ## External modules

import rebound as rb

# ## System Class
# 
# This is the most important class in the whole package.  This class allows to create the planetary system and manipulate it.

"""
#Create a simple system
#Once you create a system, a null spangler is created 

sys.set_observer(n_obs=[1,1,0],alpha_obs=0)

#Add star (by default, m = 1)
S=sys.add()

#Add planet, when an object is added, it is automatically spangled
P=sys.add("Planet",radius=0.1,m=1e-3,a=1,e=0.2)

#Add moon: orbital elements are respect to equatorial plane of the primary
M=sys.add("Planet",primary=P,radius=0.01,m=1e-7,a=0.1,e=0.01)

#Add ring system
R=sys.add("Ring",primary=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)

#If you change the number of spangles of an object the spanglers are reset
sys.update_body(R,nspangles=800)

#Each time an object is updated, the spangler should be rejoined and the simulation reset.

#Spangle 
#sys.spangle_system()

#You may check separately the properties of each object
R.spangler.plot3d()
R.spangler.plot_obs()

#Plot
sys.spangler.plot3d()
sys.spangler.plot_obs()
""";

System_doc=f"""
Creates a planetary system.

Initialization attributes:

    units: list of strings, default = ['au','msun','yr']:
        Units used in calculations following the conventions and signs of rebound.
        The order SHOULD always be MKS: length, mass, time (in that order)

Derived attributes:

    sim: Simulation:
        REBOUND Simulation object.
        
    ul, um, ut: float [SI units]
        Value of the conversion factors for each unit.
        
    G: float [ul^3/ut^2/um]
        Value of the gravitational constant.

    hashes: list:
        List of hashes of bodies in the system.
        
    stars, planets, rings: lists:
        List of the corresponding kind of object in the system.
        
    nstars, nplanets, nrings, nobservers, nbodies: int
        Number of each kind of body in the system and of all bodies in the system.

Examples:

    #Create a system
    sys=System(units=["au","msun","yr"])
    sys.sim.integrator='wahfast'
    sys.sim.dt=0.01
    
    #Add star (by default, m = 1)
    S=sys.add()

    #Add planet, when an object is added, it is automatically spangled
    P=sys.add("Planet",radius=0.1,m=1e-3,a=1,e=0.2)

    #Add moon: orbital elements are respect to equatorial plane of the primary
    M=sys.add("Planet",primary=P,radius=0.01,m=1e-7,a=0.1,e=0.01)

    #Add ring system
    R=sys.add("Ring",primary=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)

""";

class System(PrynglesCommon):
    
    def __init__(self,
                 units=['au','msun','yr'],
                ):
        
        #Rebound simulation
        self.sim=rb.Simulation()
        
        #Update rebound units
        self.update_units(units)
        
        #Bodies
        self.bodies=dict()
        self.nbodies=0
        
        #Update system
        self._update_system()
        
    def update_units(self,units):
        """Update units of the system
        """
        self.units=units
        
        #Units
        self._ul,self._um,self._ut=self.units
        self.sim.units=self.units
        
        #Canonical units of the system
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")
        self.ut=np.sqrt(self.sim.G*self.ul**3/(self.um*GSI))
        
    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        self.nparticles=len(self.sim.particles)

System.__doc__=System_doc


def add(self,kind="Star",primary=None,center="primary",**props):
    """Add an object to the system
    
    Examples:
    
        sys=System()
        S=sys.add("Star",m=2)
    
    Parameters:
    
        kind: string, default = "Star":
            Kind of object: Star, Planet, Ring (see BODY_KINDS).
    
        primary: Body, default = None:
            Primary object of the body.
            
        center: string, default = "primary":
            Center with respect to the positions are indicated.  
            Possible values: "primary", "inertial".
            When "inertial" you can provide positions of the objects using cartesian
            coordinates.

        props: dictionary:
            List of properties of the body.
            
    Returns:
        
        Body
            Body added to the system.
    """
    if kind is None:
        raise AssertionError("You must provide a valid object kind (Star, Planet, Ring).")

    if kind not in BODY_KINDS:
        raise ValueError(f"Object kind '{kind}' is not recognized.")

    exec(f"self.__body={kind}(primary=primary,**props)")
    self.bodies[self.__body.hash]=self.__body
    self.nbodies=len(self.bodies)
    
    if kind != "Ring":
        #Add body to simulation
        rb_add_options={k:v for k,v in self.__body.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}
        rb_add_options.update(hash=self.__body.hash)
        
        if primary and center=="primary":
            rb_add_options.update(primary=self.sim.particles[primary.hash])

        verbose(VERB_VERIFY,f"Adding rebound object with hash {self.__body.hash} with center {center}")
        verbose(VERB_DEEP,f"Rebound add options {rb_add_options}")
        
        self.sim.add(**rb_add_options)
        self.__body.particle=self.sim.particles[self.__body.hash]
    else:
        self.__body.particle=self.sim.particles[self.__body.primary.hash]
    
    self._update_system()
    return self.__body
    
System.add=add


def remove(self,hash):
    """Remove a body from a system.

    Parameters:
        body_hash: string
            Hash of the body to remove
    
    Notes: 
        Remove eliminate body and all the childs and the childs of the childs.

    Example:
        sys=System()
        S=sys.add(m=2)
        sys.remove(hash=S.hash)
    """
    
    if hash in self.bodies:
        verbose(VERB_SIMPLE,f"Removing object {hash} from system")

        obj=self.bodies[hash]

        #Get the list of child hashes before removing (it changes during for)
        child_hashes=list(obj.childs.keys())
        
        #Remove child objects
        for child_hash in child_hashes:
            if child_hash in self.bodies:
                self.remove(child_hash)
                
        #Remove object from simulation
        if obj.kind != "Ring":
            verbose(VERB_SIMPLE,f"Removing particle {hash} from simulation")
            self.sim.remove(hash=hash)
        
        #Remove object from childs of its primary
        if obj.primary:
            del obj.primary.childs[hash]
        
        #Remove object from bodies
        del self.bodies[hash]

        #Update system
        self._update_system()
    else:
        raise ValueError("No object with hash 'body_hash' in the system")
System.remove=remove


def spangle_system(self,n_obs=[0,0,1],alpha_obs=0):
    
    self.spanglers=[]
    for hash,body in self.bodies.items():
        print(f"Spangling body '{hash}' (kind '{body.kind}')")
        body.spangle_body()
        body.sp.set_observer(n_obs=n_obs,alpha_obs=alpha_obs)
        body.sp.set_positions(center_ecl=body.particle.xyz)
        self.spanglers+=[body.sp]
    
    #Join spanglers
    self.sp=Spangler(spanglers=self.spanglers)
    
    #Set observer
    self.sp.set_observer(n_obs=n_obs,alpha_obs=alpha_obs)
        
System.spangle_system=spangle_system

