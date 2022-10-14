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

#@external
import rebound as rb
from tqdm import tqdm
#@end

#@aliases
sci=Science
print_df=Misc.print_df
#@end

# ## System Class
# 
# This is the most important class in the whole package.  This class allows to create the planetary system and manipulate it.

#@doc:System
System_doc=f"""Creates a planetary system.

    Initialization attributes:

        units: list of strings, default = ['au','msun','yr2pi']:
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

        nbodies: int:
            Number of bodies.

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
#@end

class System(PrynglesCommon):
    
    def __init__(self,
                 filename=None,
                 units=['au','msun','yr2pi'],
                 resetable=False
                ):
        
        if filename:
            self.load_from(filename)
            return
        
        #Rebound simulation
        self.sim=rb.Simulation()
        
        #Attributes by default
        self.bodies=dict()
        self.source=None
        self.center_source=np.array([0,0,0])
        
        #Observer properties
        self.n_obs=np.array([0,0,1])
        self.alpha_obs=0  
        self.center_obs=None
        
        #Check if spangled
        self._spangled=False
        
        #Check if observer has been set
        self._observer_set=False
        self._luz_set=False
        
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
        #Check units
        if units[0] not in rb.units.lengths_SI:
            raise ValueError(f"Length unit provided '{units[0]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.lengths_SI.keys())}")
        if units[1] not in rb.units.masses_SI:
            raise ValueError(f"Mass unit provided '{units[1]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.masses_SI.keys())}")
        if units[2] not in rb.units.times_SI:
            raise ValueError(f"Time unit provided '{units[2]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.times_SI.keys())}")
        
        #Units        
        self.units=units
        self._ul,self._um,self._ut=self.units
        self.sim.units=self.units
        
        #Canonical units of the system
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")
        self.ut=np.sqrt(self.sim.G*self.ul**3/(self.um*GSI))
        
        #Update system
        self._update_system()
        
    def _get_source(self,body):
        """Get the source of light (stellar body) in the center of a body
        """
        if (body.primary is None) or (body.kind == "Star"):
            return body

        elif body.primary.kind == "Star":
            return body.primary

        else:
            return self._get_source(body.primary)

    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        self.nparticles=len(self.sim.particles)
        
    def _is_spangled(self):
        """Check if system is spangled
        """
        return True if self.sg else False
    
    def reset_state(self):
        """Reset the state of the spangler
        """
        self.sg.reset_state()
        self._observer_set=False
        self._luz_set=False

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
        print(f"System with {self.nbodies} bodies and {self.nparticles} particles (rings and disk are not particles)")
        sys.sim.status()

System.__doc__=System_doc


#@method:System
def add(self,kind="Star",primary=None,**props):
    """Add an object to the system
    
    Examples:
    
        sys=System()
        S=sys.add("Star",m=2)
    
    Parameters:
    
        kind: string, default = "Star":
            Kind of object: Star, Planet, Ring (see BODY_KINDS).
    
        primary: Body, default = None:
            Primary object of the body.
            
        props: dictionary:
            List of properties of the body.
            
    Returns:
        
        Body
            Body added to the system.
            
    Examples:
        #Add star (by default, m = 1)
        S=sys.add()

        #Add planet, when an object is added, it is automatically spangled
        P=sys.add("Planet",radius=0.1,m=1e-3,x=1,vy=0.2)

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
    props.update(dict(hash_by_kind=True))
    exec(f"self.__body={kind}(primary=primary,**props)")
    #del props["hash_by_kind"]
    
    if self.__body.bhash in self.bodies:
        raise ValueError(f"An object with hash '{self.__body.bhash}' has been already added.")
    
    self.bodies[self.__body.bhash]=self.__body
    
    if kind == "Star":
        if self.source:
            raise AssertionError("Present version of Pryngles only works with one star per system")
            
        #Add the source of light
        self.source=self.__body
    
    if kind == "Ring":
        
        #If it is a ring it does not need to be add to rebound
        self.__body.rbhash=self.__body.primary.rbhash
    
    else:
        
        #Add body to simulation
        rb_add_options={k:v for k,v in self.__body.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}
        rb_add_options.update(hash=self.__body.bhash)
        
        verbose(VERB_VERIFY,f"Adding rebound object with hash {self.__body.bhash}")
        verbose(VERB_DEEP,f"Rebound add options {rb_add_options}")
        
        #Add particle to rebound
        self.sim.add(**rb_add_options)
        
        #self.__body.particle=self.sim.particles[self.__body.hash] <- This is not convenient for pickling
        self.__body.rbhash=self.__body.bhash

    verbose(VERB_SIMPLE,f"Object '{kind}' with hash '{self.__body.bhash}' has been added.")
    #Update system
    self._update_system()
    return self.__body

#@end
System.add=add

#@test
#@end

#@method:System
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
#@end
System.remove=remove

#@test
#@end

# ## Spnagle System, Set Observer and Set Luz

#@method:System
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
    
    for bhash,body in self.bodies.items():
        
        verbose(VERB_SIMPLE,f"Spangling body '{bhash}' (kind '{body.kind}')")
        body.spangle_body()

        #Center object around its position according to rebound
        body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
        body.sg.set_positions(center_ecl=body.center_ecl)
            
        self._spanglers[bhash]=body.sg
    
        if body.kind=="Star":
            body.sg.data.source=True
            self.center_source=body.center_ecl
        
    #Join spanglers
    self.sg=Spangler(spanglers=list(self._spanglers.values()))

    #An usefule alias
    self.data=self.sg.data
    
    #Set default observer
    self.update_perspective(n_obs=self.n_obs,alpha_obs=self.alpha_obs)
    
    #Save state of the system
    if self._resetable:
        self.save_to(self._snap_file_name)
    
    #Already spangled
    self._spangled=True

def _set_observer(self,nvec=[0,0,1],alpha=0,center=None):
    """Set the position of the observer
    """
    #Only set observer if it is spangled
    if self._is_spangled():
        
        #At changing the observer, reset state
        self.sg.reset_state()
        
        #Set observer
        self.sg.set_observer(nvec=nvec,alpha=alpha,center=center)
        self.d_obs=self.sg.d_obs
        self.n_obs=self.sg.n_obs.copy()
        self.rqf_obs=self.sg.rqf_obs.copy()
        self.alpha_obs=self.sg.alpha_obs
        self.center_obs=self.sg.center_obs
    
        #Update visibility
        self.sg.update_visibility_state()
        
        #Check that observer has been set
        self._observer_set=True
        
    else:
        raise AssertionError("You must first spangle system before setting observer direction.")
        
def _set_luz_recursive(self,bhash,nluz):
    """Set light source for body and 
    """
    verbose(VERB_SIMPLE,f"Illuminating body {bhash}, with {nluz} and {self.center_source}")
    body=self.bodies[bhash]
    self.sg.set_luz(nvec=nluz,center=self.center_source,bhash=bhash)
    if body.childs:
        verbose(VERB_SIMPLE,f"Object {bhash} has childs!")
        for child_bhash in body.childs:
            self._set_luz_recursive(child_bhash,nluz)
            
def _set_luz(self):
    """Set illumination in the system.
    
    Update:
        States: illuminated, shadow, hidden_by_luz
    """
    if self._is_spangled():
        
        if not self._observer_set:
            raise AssertionError("You must first set observer before setting light.")
        
        self.bodies_illuminated=[]
        for bhash,body in self.bodies.items():
          
            if body.kind == "Star":
                verbose(VERB_SIMPLE,f"Body is a star... skipping")
                continue
                
            if body.primary.bhash in self.bodies_illuminated:
                verbose(VERB_SIMPLE,f"Primary body of {bhash}, {body.primary.bhash}, has been already illuminated")
                continue
                        
            #Get center of body
            center=body.center_ecl
                    
            #Get source and center
            verbose(VERB_SIMPLE,f"Calculating illumination for '{bhash}' coming from '{self.source.bhash}' @ {self.center_source}")            
            nluz=self.center_source-center
                    
            if body.kind == "Ring" and body.primary.kind == "Star":
                verbose(VERB_SIMPLE,f"Primary body of ring, {body.primary.bhash} is a star. All spangles will be illuminated")
                self.sg.set_luz(nvec=nluz,center=self.center_source,bhash=bhash)
                cond=(self.sg.data.bhash==bhash)
                self.sg.data.loc[cond,"unset"]=False
                self.sg.data.loc[cond,"illuminated"]=True                    
            else:                
                verbose(VERB_SIMPLE,f"Illuminating body {bhash} and all its childs")
                self._set_luz_recursive(bhash,nluz)
                self.sg.update_illumination_state()
                self.bodies_illuminated+=[bhash]
                
        self._luz_set=True
    else:
        raise AssertionError("You must first spangle system before setting light.")
        
def update_perspective(self,n_obs=None,alpha_obs=0,center_obs=None):
    """Update perspective (observer)
    """
    if n_obs is not None:
        #Update observing conditions
        self.n_obs,one=spy.unorm(n_obs)
        self.alpha_obs=alpha_obs
        self.center_obs=center_obs

    #Set observer
    self._set_observer(nvec=self.n_obs,alpha=self.alpha_obs,center=center_obs)
    self._set_luz()
    
#@end
            
System._set_luz_recursive=_set_luz_recursive
System._set_observer=_set_observer
System._set_luz=_set_luz
System.update_perspective=update_perspective
System.spangle_system=spangle_system

#@test
#@end

#@test
#@end


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


def reset(self):
    """Reset system to spangling state
    """
    if self._resetable:
        self.load_from(self._snap_file_name)
        pass
    else:
        print("System is not resetable. Use resetable = True when defining the System or when you spangle it.")
    
System.reset=reset


def integrate(self,*args,**kwargs):
    """Integrate system

    Parameters:
        *args, **kwargs:
            Mandatory (non-keyword) arguments and optional (keyword) arguments for rebound.integrate.
        
    Update:
        Integrate using integrate rebound method.
        
        Update center of each body and set positions of the spangles.
    """
    #Time of integration
    t=args[0]
    verbose(VERB_SIMPLE,"Integrating up to {t}")
    
    if self._is_spangled():
        
        #Integrate
        self.sim.integrate(*args,**kwargs)
        self.sim.move_to_com()
    
        #Update positions
        for bhash,body in self.bodies.items():
            
            #Position of the body according
            body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)

            verbose(VERB_VERIFY,f"Updating center of body {bhash} @ {body.center_ecl}")
            cond=self.sg.data.bhash==bhash
            self.sg.data.loc[cond,"center_ecl"]=pd.Series([list(body.center_ecl)]*sum(cond),dtype=object).values

        #Update positions
        self.sg.set_positions()
        
    else:
        raise AssertionError("You must first spangle system before setting positions.")
    
System.integrate=integrate


