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

# ## Aliases

sci=Science
print_df=Misc.print_df

# ## System Class
# 
# This is the most important class in the whole package.  This class allows to create the planetary system and manipulate it.

System_doc=f"""Creates a planetary system.

    Initialization attributes:

        units: list of strings, default = ['au','msun','yr']:
            Units used in calculations following the conventions and signs of rebound.
            The order SHOULD always be MKS: length, mass, time (in that order)

    Optional attributes:

        resetable: boolean, default = False:
            If True the system is resetable, namely you can reset it to the initial system.
            
        filename: string, default = None:
            File to load system.

    Derived attributes:

        sim: Class Simulation:
            Rebound Simulation object.

        ul, um, ut: float [SI units]:
            Value of the conversion factors for each unit.

        G: float [ul^3/ut^2/um]
            Value of the gravitational constant.

        bodies: dictionary:
            Bodies in the system.

        sources: dictionary:
            Bodies in the system which are sources of light.

        nbodies: int:
            Number of bodies.

        nsources: int:
            Number of sources of light.

        nparticles: int:
            Numbre of particles in rebound simulation.

        spangler: Class Spangler:
            Spangler object with all the spangles in the system.

    Examples:

        #Create a system
        sys=System(units=["au","msun","yr"])
        sys.sim.integrator='whfast'
        sys.sim.dt=0.01

        #Add star (by default, m = 1)
        S=sys.add()

        #Add planet, when an object is added, it is automatically spangled
        P=sys.add("Planet",radius=0.1,m=1e-3,a=1,e=0.2)

        #Add moon: orbital elements are respect to ecliptic system
        M=sys.add("Planet",primary=P,radius=0.01,m=1e-7,a=0.1,e=0.01)

        #Add ring system
        R=sys.add("Ring",primary=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)

""";

class System(PrynglesCommon):
    
    def __init__(self,
                 filename=None,
                 units=['au','msun','yr'],
                 resetable=False
                ):
        
        if filename:
            self.load_from(filename)
            return
        
        #Rebound simulation
        self.sim=rb.Simulation()
        
        #Attributes by default
        self.bodies=dict()
        self.sources=dict()
        
        #Observer properties
        self.n_obs=[0,0,1]
        self.alpha_obs=0  
        
        #Check if spangled
        self._spangled=False

        #Initialize spangler object
        self.sg=None
        
        #Is the system resetable?
        self._resetable=resetable
        if self._resetable:
            #Create temporary file
            self._snap_file_name = "/tmp/pryngles-system.pkl"
        
        #Update rebound units
        self.update_units(units)
        
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
        
        #Update system
        self._update_system()
        
    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        self.nsources=len(self.sources)
        self.nparticles=len(self.sim.particles)
        
    def _is_spangled(self):
        return True if self.sg else False
        
    def save_to(self,filename):
        """Save system from file
        
        Parameters:
            filename: string:
                Path to file where the object will be pickled.
                
        Result:
            File 'filename' for regular object and 'filename.rbin' for rebound simulation
        """
        #Rebound file
        rb_filename=filename+".rbin"

        #Save rebound state
        verbose(VERB_SIMPLE,"Saving rebound simulation")
        self.sim.save(rb_filename)

        #Since rebound have ctypes it cannot be pickled
        del self.sim

        #Pickle system
        PrynglesCommon.save_to(self,filename)

        #Load again rebound
        self.sim=rb.Simulation(rb_filename)

    def load_from(self,filename):
        """Load system from filename
                
        Parameters:
            filename: string:
                Path to file where the object will be pickled.
                There to be 2 files: 'filename' (with the regular object) and filename.rbin with 
                rebound simulation.
        """
        #Rebound file
        rb_filename=filename+".rbin"

        #Load system
        self=PrynglesCommon.load_from(self,filename)

        #Load rebound
        verbose(VERB_SIMPLE,"Loading rebound simulation")
        self.sim=rb.Simulation(rb_filename)
        
    def status(self):
        print(f"System with {self.nbodies} bodies, {self.nsources} sources and {self.nparticles} particles")
        sys.sim.status()

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
            
    Examples:
        #Add star (by default, m = 1)
        S=sys.add()

        #Add planet, when an object is added, it is automatically spangled
        P=sys.add("Planet",radius=0.1,m=1e-3,a=1,e=0.2)

        #Add moon: orbital elements are respect to ecliptic system
        M=sys.add("Planet",primary=P,radius=0.01,m=1e-7,a=0.1,e=0.01)

        #Add ring system
        R=sys.add("Ring",primary=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)        
        
    """
    if kind is None:
        raise AssertionError("You must provide a valid object kind (Star, Planet, Ring).")

    if kind not in BODY_KINDS:
        raise ValueError(f"Object kind '{kind}' is not recognized.")

    #Create body
    exec(f"self.__body={kind}(primary=primary,**props)")
    
    if self.__body.bhash in self.bodies:
        raise ValueError(f"An object with hash '{self.__body.bhash}' has been already added.")
    
    self.bodies[self.__body.bhash]=self.__body
    
    if kind == "Star":
        #Add a source
        self.sources[self.__body.bhash]=self.__body
    
    if kind == "Ring":
        
        #If it is a ring it does not need to be add to rebound
        self.__body.rbhash=self.__body.primary.rbhash
    
    else:
        
        #Add body to simulation
        rb_add_options={k:v for k,v in self.__body.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}
        rb_add_options.update(hash=self.__body.bhash)
        
        #Check if positions are given with respect to primary
        if primary and center=="primary":
            rb_add_options.update(primary=self.sim.particles[primary.bhash])

        verbose(VERB_VERIFY,f"Adding rebound object with hash {self.__body.bhash} with center {center}")
        verbose(VERB_DEEP,f"Rebound add options {rb_add_options}")
        
        #Add particle to rebound
        self.sim.add(**rb_add_options)
        
        #self.__body.particle=self.sim.particles[self.__body.hash] <- This is not convenient for pickling
        self.__body.rbhash=self.__body.bhash

    verbose(VERB_SIMPLE,f"Object '{kind}' with hash '{self.__body.bhash}' has been added.")
    #Update system
    self._update_system()
    return self.__body
    
System.add=add


def remove(self,bhash):
    """Remove a body from a system.

    Parameters:
        bhash: string
            Hash of the body to remove
    
    Notes: 
        Remove eliminate body and all the childs and the childs of the childs.

    Example:
        sys=System()
        S=sys.add(m=2)
        sys.remove(bhash=S.bhash)
    """
    
    if bhash in self.bodies:
        verbose(VERB_SIMPLE,f"Removing object {bhash} from system")

        obj=self.bodies[bhash]

        #Get the list of child hashes before removing (it changes during for)
        child_hashes=list(obj.childs.keys())
        
        #Remove child objects
        for child_hash in child_hashes:
            if child_hash in self.bodies:
                self.remove(child_hash)
                
        #Remove object from simulation
        if obj.kind != "Ring":
            verbose(VERB_SIMPLE,f"Removing particle {bhash} from simulation")
            self.sim.remove(hash=bhash)
        
        #Remove object from childs of its primary
        if obj.primary:
            del obj.primary.childs[bhash]
        
        #Remove object from bodies
        del self.bodies[bhash]

        #Update system
        self._update_system()
    else:
        raise ValueError(f"No object with hash '{bhash}' in the system")

System.remove=remove


def integrate(self,*args,**kwargs):
    """Integrate system
    
    *args, **kwargs:
        Mandatory (non-keyword) arguments and optional (keyword) arguments for rebound.integrate.
    """
    self.sim.integrate(*args,**kwargs)
    
System.integrate=integrate


# ## Spangle System

def spangle_system(self):
    """Generate the spangles of the objects in the system
    
    Attributes created:
        
        spanglers: dictionary of Spangler objects:
            Spangler corresponding to each object in the system.
            
        sp: Spangler:
            Spangler corresponding to all system.
            
    Result:
        
        This method create the spangler of the system

    """
    
    self._spanglers=dict()
    source=1
    for bhash,body in self.bodies.items():
        
        verbose(VERB_SIMPLE,f"Spangling body '{bhash}' (kind '{body.kind}')")
        body.spangle_body()

        if body.kind=="Star":
            body.sg.data.source=source
            source+=1
        
        #Center object around its position according to rebound
        body.sg.set_positions(center_ecl=self.sim.particles[body.rbhash].xyz)
            
        self._spanglers[bhash]=body.sg
    
    #Join spanglers
    self.sg=Spangler(spanglers=list(self._spanglers.values()))
    
    #Add column for controlling information on sources
    for source in range(1,self.nsources+1):
        for source_state in SPANGLER_SOURCE_STATES.keys():
            self.sg.data[source_state+f"_{source}"]=0
    
    #Save state of the system
    if self._resetable:
        self.save_to(self._snap_file_name)
    
    #Already spangled
    self._spangled=True

System.spangle_system=spangle_system


# ### Miscelaneous methods

def update_body(self,body,**props):
    """Update properties of a body in the system
    
    Parameters:
        body: string or Body:
            Body to update
        
        props: dict:
            Dictionary with properties of the object
    """
    #Update spangling?
    if self._is_spangled():
        raise AssertionError("After spangling you cannot update the properties of the bodies.  Please rebuild the system")

    #Update body properties
    if isinstance(body,Body):
        body.update_body(**props)
    elif body in self.bodies:
        body=self.bodies[body]
        lkind=body.kind.lower()
        exec(f"body.update_{lkind}()")
    else:
        raise AssertionError("You are trying to update a body ({body}) which is not in the system")
        
    #Check if among props there is any property related to position
    if any(k in props for k in REBOUND_ORBITAL_PROPERTIES):
        raise ValueError(f"You cannot update an orbital property {props} without compromising the full simulation. Rebuild the system from scratch.")
        
System.update_body=update_body


def set_observer(self,nvec=[0,0,1],alpha=0):
    """Set the position of the observer
    """
    #Only set observer if it is spangled
    self.n_obs=nvec
    self.alpha_obs=alpha
    if self._is_spangled():
        #Set observer
        self.sg.set_observer(nvec=self.n_obs,alpha=self.alpha_obs)
    else:
        raise AssertionError("You must first spangle system before setting observer direction.")
        
System.set_observer=set_observer


def reset(self):
    """Reset system to spangling state
    """
    if self._resetable:
        self.load_from(self._snap_file_name)
        pass
    else:
        print("System is not resetable. Use resetable = True when defining the System or when you spangle it.")
    
System.reset=reset


# ## Illumination

def update_illumination(self):
    """Determine the visibility conditions of the spangles
    
    Update n_luz differently according to body hash
    """
    #Calculate the hull of the bodies in the system
    for body in self.bodies.values():
        #Body information
        bhash=body.hash
        bkind=body.kind
    
System.update_illumination=update_illumination

