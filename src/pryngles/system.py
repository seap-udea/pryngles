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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
from pryngles import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External required packages
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import rebound as rb
from tqdm import tqdm
import spiceypy as spy
import numpy as np
from copy import deepcopy
from anytree import NodeMixin,RenderTree

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BODY_KINDS=[]

"""
These are the default attributes for any body.
"""
BODY_DEFAULTS=dict()
BODY_DEFAULTS.update(odict(
    
    name=None,
    name_by_kind=False,
    source=None,
    
    #Orbit
    #No default orbital properties required
    m=1,

    #Physics
    radius=1,
    prot=1,
    i=0, #Inclination of the rotational axis
    roll=0,
    alpha=0, #Zero meridian
    q0=0,
    
    #Optics
    nspangles=1000,
    spangle_type=SPANGLE_SOLID_ROCK,
    shape="sphere",
    geometry_args=dict(),
    seed=0,
    preset=True,
    spangles=None,
    normarea=1,

    #Basic optical properties
    albedo_gray_normal=1,
    tau_gray_optical=0,
    
    #Legacy
    primary=None,
    optics=dict(),
    orbit=dict(),
    physics=dict(),
))
BODY_KINDS=[]

"""
These are the default attributes for bodies of the kind 'Star'.
"""
STAR_DEFAULTS=deepcopy(BODY_DEFAULTS)
STAR_DEFAULTS.update(odict(

    #Orbit: update
    #No default orbital properties required
    
    #Physics: update
    #Same as Body
    
    #Optical properties: update
    limb_coeffs=[],
    spangle_type=SPANGLE_STELLAR,
    shape="sphere",
))
BODY_KINDS+=["Star"]

"""
These are the default attributes for bodies of the kind 'Planet'.
"""
PLANET_DEFAULTS=deepcopy(BODY_DEFAULTS)
PLANET_DEFAULTS.update(odict(
    
    #Orbit: update
    #No default orbital properties required
    
    #Physics: update
    #Same as Body
    radius=0.1,
    
    #Optical: update
    spangle_type=SPANGLE_SOLID_ROCK,
    geometry="sphere",
))
BODY_KINDS+=["Planet"]

"""
These are the default attributes for bodies of the kind 'Ring'.
"""
RING_DEFAULTS=deepcopy(BODY_DEFAULTS)
RING_DEFAULTS.update(odict(

    #Orbit: update
    #Same as Body altough ring has not orbit properties
    
    #Physics: update
    #Same as Body
    fi=1.5,
    fe=2.0,
    taur=0.4,
    
    #Optics: update
    spangle_type=SPANGLE_GRANULAR,
    shape="ring",
))
BODY_KINDS+=["Ring"]

"""
These are the default attributes for bodies of the kind 'Observer'.
"""
OBSERVER_DEFAULTS=deepcopy(BODY_DEFAULTS)
OBSERVER_DEFAULTS.update(odict(
    lamb=0,
    beta=0,
))
BODY_KINDS+=["Observer"]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Body
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Body(Orbody):
    """A general body.  This class is not intended to be used independently, just for inheritance purposes.
        
    Initialization attributes:
    
        kind : string:
            One of the kind of bodies defined in the package (see _BODY_KINDS)
            Defined objects are: "Star", "Planet", "Ring".
    
        defaults : OrderedDict:
            Dictionary with the properties of the object.
    
        parent: Class Body:
            Object in the center of the orbit of this body.
    
        **properties: dicitionary:
            Specification of the body properties.  All objects of the class Body has the following
            properties by default:
            
            name: string, default = None:
                Name of the object, ie. a unique string identifying the object.  It can be provided
                by the user or automatically set by the initializer using a unique hash 
                (see hash Python function).
    
            orbital properties: 
                Object with the orbital properties of the body (eg. orbit.m is the mass)
                see each specific Body definition for attributes.
                orbit must be compatible with rebound.
    
                    m: float [rebound mass units], default = 1:
                        Mass of the body.  If m = 0 the body does not produce gravitation.
    
            physical properties:
    
                Object with the physical properties of the body (eg. physics.radius)
                see each specific Body definition for attributes.
    
                    radius: float [rebound length units], default = 1:
                        Radius of the body.
    
                    prot: float [rebound time units], default = 1:
                        Period of rotation of the body.
    
                    i: float [rad], default = 0:
                        Inclination of the body equator with respect to the ecliptic plane.
    
                    roll: float [rad], default = 0:
                        Roll angle.  This is the angle with respect to ecliptic x-axis in which 
                        the normal to the object equatorial plane is rotated.
    
                    alpha_equ: float [rad], default = 0:
                        Longitude of the zero meridian of the object.
    
                    q0: float [ut], default = 0:
                        Initial longitude for zero meridian.
    
            optical properties:
    
                Object with the optical properties of the body (eg. physics.lamb_albedo)
                see each specific Body definition for attributes.
    
                    nspangles: int, default = 1000:
                        Number of spangles on which the object will be discretized.
                        
                    spangle_type: int, default = SOLID_SPANGLE:
                        Type of spangles of the body.
                        
                    preset: boolean, default = True:
                        If True spangle object from a preset.
    
    Derived attributes:
    
            wrot: float [rad/rebound time units]:
                Rotational angular velocity.
    
            n_equ: array(3):
                Rotational axis vector in the ecliptic system.
        
    Secondary attributes:
    
        childs: list
            List with child bodies (bodies which have this body) as the center.
    
    Public methods:
    
        update_body(**props):
            Update a given set of properties.
            
    Examples:
    
        Create a body with None parent and name = 'B':
        
            B=Body("Body",BODY_DEFAULTS,None,name='B',m=2,c=2)
            
        Create a body having parent the Body "B" defined before:
             
            C=Body("Body",BODY_DEFAULTS,B,name="C")
    
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    def __init__(self,kind,defaults,parent,**props):

        #Kind, parent and child attributes
        self.kind=kind
        self.__defaults=defaults
        
        #Prepare key attributes
        self.sg=None
        self._spangles = None
        self._normarea = 1

        #Name of the object
        if 'name' in props:
            name=self.name=str(props["name"])
        elif 'name_by_kind' in props:
            name=self.name=self.kind
        else:
            name=self.name=str(hash(self))

        #Legacy
        if 'primary' in props:
            parent=props["primary"]
        if 'optics' in props:
            props.update(props["optics"])
        if 'orbit' in props:
            props.update(props["orbit"])
        if 'physics' in props:
            props.update(props["physics"])

        #Update childs and parent
        if parent is not None:
            if not isinstance(parent,Body):
                raise AssertionError(f"Parent is not a valid Object: {type(parent)}, {isinstance(parent,Body)}")
            else:
                self.parent=parent
                parent._update_childs(self)

        #Update parent and childs        
        self._update_parent(parent)
        self._update_childs()

        #Update default properties
        self.__dict__.update(defaults)
        #Set name
        self.name=name
        #Update body
        self.update_body(**props)
    
    def update_body(self,**props):
        """Update properties of the Body.
        
        Parametes:
            **props: dictionary:
                Properties to update. The current object is updated with new 
                values provided in this new object
                
        Example:
            B.update_body(m=2)
                This only update the attribute m of orbit.
        """
        for prop in props:
            if prop in self.__defaults or prop in REBOUND_ORBITAL_PROPERTIES:
                self.__dict__[prop]=props[prop]
            else:
                raise ValueError(f"Property {prop} not identified in object {self.kind}")
                
        self.elements={k:v for k,v in self.__dict__.items() if k in REBOUND_ORBITAL_PROPERTIES}

        """We introduce this code to allow the initialization of a
        body with the coordinates of the spangles instead of
        initializing them
        """
        try:
            nspangles = len(self.nspangles["spangles"])
            self._spangles = np.array(self.nspangles["spangles"])
            self._normarea = self.nspangles["normarea"] if 'normarea' in self.nspangles.keys() else 1
            self.nspangles = nspangles
        except:
            pass
        verbose(VERB_VERIFY,"Updating Body")
        self._update_properties()
    
    def _update_childs(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=dict()
        if child is not None:
            verbose(VERB_VERIFY,f"Add child {child.name} to body {self.kind} ({self.name})")
            self.childs[child.name]=child
            
    def _update_parent(self,parent=None):
        if 'parent' not in self.__dict__:
            if parent:
                verbose(VERB_VERIFY,f"Add parent {parent.name} to body {self.kind} ({self.name})")
            self.parent=parent
        elif parent is not None:
            verbose(VERB_VERIFY,f"Add parent {parent.name} to body {self.kind} ({self.name})")
            self.parent=parent
            parent._update_childs(self)
    
    def _update_properties(self):
        verbose(VERB_VERIFY,"Updating properties of Body")
        #Rotational angular velocity
        self.wrot=2*np.pi/self.prot
        #Rotation axis
        self.n_equ=sci.cartesian([1,self.roll,90*Consts.deg-self.i])
    
    def show_tree(self):
        print(RenderTree(self))
        

    def spangle_body(self):

        """
        Spangle the surface of the body
        """
        
        #Create spangler
        self.sg=Spangler(
            nspangles=self.nspangles,
            spangles=self._spangles,
            normarea=self._normarea,
            name=self.name,
            n_equ=self.n_equ,
            alpha_equ=self.alpha,
            w=self.wrot,
            q0=self.q0,
        )
        
        #Populate spangler
        self.sg.populate_spangler(
            shape=self.shape,
            spangle_type=self.spangle_type,
            scale=self.radius,
            seed=self.seed,
            preset=self.preset,
            **self.geometry_args,
        )
        
        #Additional properties in the Spangler DataFrame
        if self.kind=="Star":
            self.sg.data.source=True
        
        self.sg.set_observer()
        self.sg.set_luz()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Star
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Star(Body):
    """A star.

    Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.

            If None the object is the center of the orbit specification for other objects.
            
            Object parent for a star should be another star.
        
        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:
            
                limb_coeffs: list [adimensional], default = []:
                    List of limb darkening fit coefficients.  See Science.calc_limbdarkening.

                    Models in: https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html
                    Coefficients available at: https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt
                    
                spangle_type: int, default = STAR_SPANGLE:
                    Type of spangles

    Derived attributes:
    
    Methods:
    
        update_body(**pars):

            This method compute some derived attributes like.

    Notes:

        See Body class documentation.
    
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Star",STAR_DEFAULTS,parent,**props)

        #Check parent
        if self.parent is not None:
            if self.parent.kind!="Star":
                raise ValueError(f"Only another Star can be the parent of a Star (you provided {self.parent.kind})")

        self._update_star_properties()
        
    def _update_star_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            norm_limb_darkening: float:
                Limb darkening function normalization.
                Requires: limb_coefs.

        """
        verbose(VERB_VERIFY,"Updating properties of Star")
        
        #Compute limbdarkening at r = 0 to initialize normalization constant
        sci.limb_darkening(0,self.limb_coeffs)
        
        #Store limb darkening normalization
        self.norm_limb_darkening=SCIENCE_LIMB_NORMALIZATIONS[hash(tuple(self.limb_coeffs))]
        
    def update_star(self,**props):
        """General update propeties of the Star
        """
        verbose(VERB_VERIFY,"Updating star")
        
        Body.update_body(self,**props)
        self._update_star_properties()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Planet
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Planet(Body):
    """A planet.

    Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.

        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:
            
                x,y,z: float [ul], default = 1.0, 0.0, 0.0:
                    Initial position of the body.

                vy: float, default = 0.0, 1.0, 0.0:
                    Intitial velocity of the body
        
    Derived attributes:
        None.
    
    Notes:

        See Body class documentation.
    
    """
    
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Planet",PLANET_DEFAULTS,parent,**props)
        
        #Check parent
        if self.parent is None:
            raise ValueError(f"Parent not provided and it is mandatory for {self.kind}.")
        
        #Update properties
        self.update_planet(**props)

    def _update_planet_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            norm_limb_darkening: float:
                Limb darkening function normalization.
                Requires: limb_coefs.

        """
        verbose(VERB_VERIFY,"Updating Planet properties")
        
    def update_planet(self,**pars):
        verbose(VERB_VERIFY,"Updating Planet")
        Body.update_body(self,**pars)
        self._update_planet_properties()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Ring
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Ring(Body):
    """Class Ring.
    
Initialization attributes:
        
        parent: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.

        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:

            fi: float [adimensional], default = 1:
                Fraction of the radius of the parent object where ring stars.

            fe: float [adimensional], default = 1:
                Fraction of the radius of the parent object where ring ends.

            albedo_gray_normal: float. default = 1: 
                Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.

            tau_gray_optical: float. default = 0:
                Gray (wavelength indpendent) Optical depth of the spangle.  
                If 0 the spangle is entirely opaque to all wavelength, despite its type.            

    Derived attributes:
    
        ri: float:
            Radius of the inner border of the ring in units of the parent radius.

        re: float:
            Radius of the outer border of the ring in units of the parent radius.
            
    Notes:

        See Body class documentation.
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Ring",RING_DEFAULTS,parent,**props)
        
        #Check parent
        if self.parent is None:
            raise ValueError(f"Parent not provided and mandatory for {self.kind}.")
        
        #Update properties
        self.update_ring(**props)

    def _update_ring_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            ri, re: float:
                Radius of the inner (outer) border of the ring in units of the parent radius.
                Requires: limb_coefs.
                
            radius: float:
                Object radius.
                
            geometry_args: dictionary:
                
        """
        verbose(VERB_VERIFY,"Updating Ring properties")
    
        #Update radius
        self.ri=self.fi*self.parent.radius
        self.re=self.fe*self.parent.radius
        self.radius=self.re
        
        #Update geometry args for spangling purposes
        self.geometry_args=dict(ri=self.ri/self.re)
        
    def update_ring(self,**pars):
        verbose(VERB_VERIFY,"Updating Ring")
        Body.update_body(self,**pars)
        self._update_ring_properties()   



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Observer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Observer(Body):
    """This class is intended only for legacy purposes.
    """
    def __init__(self,
                 parent=None,
                 **props
                ):
        Body.__init__(self,"Observer",OBSERVER_DEFAULTS,parent,**props)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class System
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class System(PrynglesCommon):
    """Creates a planetary system.
    
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
    
            #Add moon: orbital elements are respect to central object and w.r.t. to ecliptic system
            M=sys.add("Planet",parent=P,radius=0.01,m=1e-7,a=0.1,e=0.01)
    
            #Add ring system
            R=sys.add("Ring",parent=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)
    
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def __init__(self,
                 filename=None,
                 units=['au','msun','yr2pi'],
                 resetable=False
                ):
        
        if filename:
            self.load_from(filename)
            return
        
        #Rebound simulation
        self.sim=None
        self._simulated=False
        
        #Attributes by default
        
        #List of bodies in the system
        self.bodies=odict()
        
        #Root of the tree of bodies
        self.root=None
        
        #Center of the light-source in the system
        self.source=None
        self.center_root=np.array([0,0,0])
        
        #Orbital configuration
        self.orbital_configuration=None
        
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
        
        #By default spangle scatterers
        """
        This is the list of the class of scatterers used to calculate the scattering in
        different types of spangles.
        
        The structure of this dictionary is:
        
            key: integer (or enumerator):
                This is the column spangle_type in the spangler DataFrame.
                
            item: tuple (2):
                Component 1: 
                    Class of scatterer.
                Component 2: 
                    Dictionary mapping the initialization properties of the scatterer to columns
                    in the spangler DataFrame.
        
        Example of item:
        
            SPANGLE_ATMOSPHERIC:(GrayAtmosphere,dict(AS=1.0))
            
                This means that for spangles of the type SPANGLE_ATMOSPHERIC Pryngles will 
                instantiate an object of the class LambertianGrayAtmosphere.  This class have a
                single parameter, the spherical albedo AS.  The dictionary means that when 
                instantiating the object the column "albedo_gray_spherical" will be used to 
                initialize the object.                
        """
        self.spangle_scatterers={
            SPANGLE_ATMOSPHERIC:(GrayAtmosphere,dict(AS=1.0)),
            SPANGLE_GRANULAR:(GraySurface,dict(AL=1.0)),
            SPANGLE_LIQUID:(GraySurface,dict(AL=1.0)),
            SPANGLE_SOLID_ICE:(GraySurface,dict(AL=1.0)),
            SPANGLE_SOLID_ROCK:(GraySurface,dict(AL=1.0)),
            SPANGLE_GASEOUS:(BlackSurface,dict()),
            SPANGLE_STELLAR:(BlackSurface,dict()),
        }
        
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
        #self.sim.units=self.units
        
        #Canonical units of the system
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")

        #Compute the units of time
        sim=rb.Simulation()
        sim.units=self.units
        self.G=sim.G
        self.ut=np.sqrt(self.G*self.ul**3/(self.um*GSI))
        
        #Update system
        self._update_system()
        
    def _get_source(self,body):
        """Get the source of light (stellar body) in the center of a body
        """
        if (body.parent is None) or (body.kind == "Star"):
            return body

        elif body.parent.kind == "Star":
            return body.parent

        else:
            return self._get_source(body.parent)

    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        if self._simulated:
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
        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Save rebound state
            verbose(VERB_SIMPLE,"Saving rebound simulation")
            self.sim.save(rb_filename)

            #Since rebound have ctypes it cannot be pickled
            del self.sim
            self._simulated=True

        #Pickle system
        PrynglesCommon.save_to(self,filename)

        if self._simulated:
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
        #Load system
        self=PrynglesCommon.load_from(self,filename)

        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Load rebound
            verbose(VERB_SIMPLE,"Loading rebound simulation")
            self.sim=rb.Simulation(rb_filename)
        
    def status(self):
        if self._simulated:
            print(f"System with {self.nbodies} bodies and {self.nparticles} particles (rings and disk are not particles)")
            self.sim.status()
        else:
            print(f"Simulation for this system has not been yet initialized. Use System.initialize_simulation()")


    def update_scatterers(self, force=False, verbosity=VERB_NONE):
        """Update the scatterers of the spangles. For each spangle
        type it took the spangle_scatterers list specification and
        reassign it to all spangles in simulation.

        Optional parameters:
        
           force: boolean, default = True:
               If False it update only the unassigned spangles (this
               is to avoid initialization)
        

        """
        if not self._spangled:
            raise AssertionError("You need to spangle the system before updating the scatterers.")

        Misc.elapsed_time(show=False,verbosity=verbosity)        
        
        #Update scatterer only for the non-assigned one
        cond=(~self.data.hidden)
        cond=(cond) if force else (cond)&(self.data.scatterer=="")

        verbose(verbosity,f"Updating scatterer for {cond.sum()} spangles")
        spangle_types=dict()
        for index in self.data[cond].index:
        
            #Get spangle
            spangle=self.data.loc[index]

            #Get spangle sype
            spangle_type=spangle["spangle_type"]
                
                
            #Get scatterer class and options description
            spangle_scatterer,scatterer_parameters=self.spangle_scatterers[spangle_type]
            
            if spangle_type not in spangle_types.keys():
                verbose(verbosity,f"Updating scatterer for spangle type: '{spangle_type}' with scatterer '{spangle_scatterer}' and parameters {scatterer_parameters}")
                spangle_types[spangle_type]=True

            #Instantiate object of scatterer and save hash into DataFrame
            self.data.loc[index,"scatterer"]=spangle_scatterer(**scatterer_parameters).hash
            self.data.loc[index,"scatterer_parameters"]=scatterer_parameters,

        Misc.elapsed_time(show=True,msg="Update scatterers time",verbosity=verbosity)
            
    def add(self,kind="Star",parent=None,**props):
        """Add an object to the system
        
        Examples:
        
            sys=System()
            S=sys.add("Star",m=2)
        
        Parameters:
        
            kind: string, default = "Star":
                Kind of object: Star, Planet, Ring (see BODY_KINDS).
        
            parent: Body, default = None:
                Parent object of the body.
                
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
            M=sys.add("Planet",parent=P,radius=0.01,m=1e-7,a=0.1,e=0.01)
    
            #Add ring system
            R=sys.add("Ring",parent=P,fi=1.5,fe=2.5,albedo_gray_normal=0.5,tau_gray_optical=3)        
            
        """
        if kind is None:
            raise AssertionError("You must provide a valid object kind (Star, Planet, Ring).")
    
        if kind not in BODY_KINDS:
            raise ValueError(f"Object kind '{kind}' is not recognized.")
    
        #Legacy
        if 'primary' in props:
            parent=props["primary"]
        if kind=="Observer":
            parent=self.root
        
        #Default parameters
        if self.root:
            if (kind!="Star") and (parent is None):
                parent=self.root

            """
            if kind=="Planet":
                if "m" not in props:
                    props["m"]=0.1*parent.m
                if "radius" not in props:
                    props["radius"]=0.5*parent.radius
                if "a" not in props:
                    if "a" in parent.__dict__:
                        props["a"]=0.5*parent.a
            if kind=="Planet":
                for prop in ['m','radius']:
                    if prop not in props:
                        raise AssertionError(f"No value for '{prop}' was provided")
            """
                    
        #Create body
        props.update(dict(name_by_kind=True))
        self.__body=eval(f"{kind}(parent=parent,**props)")
        
        if self.__body.name in self.bodies:
            raise ValueError(f"An object with name '{self.__body.name}' has been already added.")
        
        #If we have a root object and no parent has been provided
        if not parent:
            if self.root:
                raise ValueError(f"A root object alread exist in the system ({self.root.name}) and you do not provided a parent body for {self.__body.name}.")
            else:
                self.root=self.__body
                verbose(VERB_SIMPLE,f"Setting the root object as {self.root.name}")
            
        self.bodies[self.__body.name]=self.__body
        
        #Create the shined body tree
        self.__body.shined=[]
    
        #Update system
        self._update_system()
        
        #Set the source of the object
        if self.__body.source:
            #Check that the source is a body
            if not isinstance(self.__body.source,Body):
                raise ValueError(f"The source of body must be an actual Body.")        
            #Check that the source is among the bodies
            if self.__body.source.name not in self.bodies:
                raise ValueError(f"The source of body {self.__body.name} is not among system bodies {list(self.bodies.keys())}.")
            #Check that the source be a star
            if self.__body.source.kind!="Star":
                raise ValueError(f"The source of body {self.__body.name} must be a Star.  You set {self.__body.source.name} which is a {self.__body.source.kind}.")
            self.__body.source.shined+=[self.__body.name]
        else:
            if self.__body.kind=="Star":
                self.__body.source=self.__body
            elif self.__body.parent==self.root:
                self.__body.source=self.root
            else:
                self.__body.source=self.__body.parent.source
        self.__body.source.shined+=[self.__body.name]
        
        verbose(VERB_SIMPLE,f"Object '{kind}' with name '{self.__body.name}' has been added.")
        return self.__body
    
    
    def initialize_simulation(self,orbital_tree=None, test=False, hinb=True, **rebound_options):
        """Initialize rebound simulation using a given orbital tree.
        
        Parameters:
            orbital_tree: list of pairs, default = None:
                A scheme showing how the bodies in the system are organized as in a 
                hierarchical N-body system (see OrbitUtil.build_system).
                
                Examples:
                    Simple system: star (S), planet (P):
                        orbital_tree = [S,P]
                    
                    System with two planets: star (S), planet 1 (P1), planet 2 (P2):
                        orbital_tree = [[S,P1],P2]
                        
                    System with moon: star (S), planet (P), moon (M):
                        orbital_tree = [S,[P,M]]
                        
                    System with two planets and moons: star (S), planet 1 (P1), moon planet 1 (M), planet 2 (P2):
                        orbital_tree = [[S,[P1,M]],P2]
                        
                    System with two stars and one planet per star:
                        orbital_tree = [[S1,PS1],[S1,PS2]]

                 When orbital_tree = None the tree will be created using the list of bodies.

            hinb: boolean, default = True:
                 If Ture Initialize the system as a hierarchical n-body system.
                 Otherwise initialize the simulation with the orbital elements provided.
        
            rebound_options: dictionary:
                 All the additional options you want to add to rebound add method.
        
        Return:
            orbit: object Orbit:
                Object containing the hierarchical N-body system.
        
        Update:
            self.sim: Rebound Simulation:
                Simulation of the system.
                
        """
        
        #Compile orbital configuration
        if orbital_tree is None:
            i=0
            for name,body in odict(reversed(list(self.bodies.items()))).items():
                if body == self.root:
                    continue
                if body.kind == "Ring":
                    continue
                if i == 0:
                    self.orbital_tree=body
                else:
                    self.orbital_tree=[body,self.orbital_tree]
                i+=1
            self.orbital_tree=[self.root,self.orbital_tree]
        else:
            self.orbital_tree=orbital_tree

        self.n_bodies_tree=len(list(Misc.flatten(self.orbital_tree)))
        
        #Set the rebound hash of all bodies
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                body.rbhash=body.parent.name
            else:
                body.rbhash=body.name
        
        #Check that all bodies in system is in the orbital tree
        bodies=list(Misc.flatten(self.orbital_tree))
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                continue
            if body not in bodies:
                raise AssertionError(f"Body '{name}' is in System but not in orbital tree.")

        #Code to test initialization
        if test:
            def show_tree(tree):
                if not isinstance(tree,list):
                    body=tree
                    print(f"Body: {body.name}")
                    print(f"Orbital elements: {body.elements}")
                else:
                    for body in tree:
                        show_tree(body)
            print(f"An orbital tree was build with {self.n_bodies_tree} bodies")
            print(f"Tree: {self.orbital_tree}")
            show_tree(self.orbital_tree)
            return self.orbital_tree
            
        if hinb:
            #Build hierarchical N-body system
            verbose(VERB_VERIFY,f"Simulation: building hierarchical n-body system with tree {self.orbital_tree}")
            orbit,pelements=OrbitUtil.build_system(self.orbital_tree,self.units)
            orbit.calculate_orbit()
            orbit.sim.move_to_com()
        else:
            #Add objects to simulation as they were provided in the adding procedures
            verbose(VERB_VERIFY,f"Simulation: adding bodies as they were provied")
            orbit=Orbit(units=self.units)
            for body in bodies:
                verbose(VERB_VERIFY,f"Adding body: {body.name}, elements: {body.elements}")
                orbit.sim.add(**body.elements)
            
        #Initialize simulation
        self.sim=rb.Simulation(**rebound_options)
        self.sim.units=self.units
        
        #Add particles to simulation
        self.sim.hash_name=dict()
        for i,p in enumerate(orbit.sim.particles):
            self.sim.add(
                hash=bodies[i].name,
                m=bodies[i].m,
                x=p.x,y=p.y,z=p.z,
                vx=p.vx,vy=p.vy,vz=p.vz
            )
            self.sim.hash_name[str(self.sim.particles[bodies[i].name].hash)]=bodies[i].name
        self.sim.orbit=orbit
        self._simulated=True
        self._update_system()
        
        return orbit
    
    def remove(self,name):
        """Remove a body from a system.
    
        Parameters:
            name: string
                Hash of the body to remove
        
        Notes: 
            Remove eliminate body and all the childs and the childs of the childs.
    
        Example:
            sys=System()
            S=sys.add(m=2)
            sys.remove(name=S.name)
        """
        
        if name in self.bodies:
            verbose(VERB_SIMPLE,f"Removing object {name} from system")
    
            obj=self.bodies[name]
    
            #Get the list of child hashes before removing (it changes during for)
            child_hashes=list(obj.childs.keys())
            
            #Remove child objects
            for child_hash in child_hashes:
                if child_hash in self.bodies:
                    self.remove(child_hash)
                    
            #Remove object from Rebound simulation
            if obj.kind != "Ring":
                if self._simulated:
                    if self.nparticles:
                        verbose(VERB_SIMPLE,f"Removing particle {name} from simulation")
                        self.sim.remove(hash=name)
            
            #Remove object from childs of its parent
            if obj.parent:
                del obj.parent.childs[name]
            
            #Remove object from bodies
            del self.bodies[name]
    
            #Update system
            self._update_system()
        else:
            raise ValueError(f"No object with hash '{name}' in the system")
    
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
        if not self._simulated:
            raise AssertionError("Before spangling the system you must initialize the simulation: System.initialize_simulation().")

        Misc.elapsed_time(show=False,verbosity=VERB_NONE)
        self._spanglers=dict()
        
        #Add spangles
        for name,body in self.bodies.items():
            
            verbose(VERB_SIMPLE,f"Spangling body '{name}' (kind '{body.kind}')")
            body.spangle_body()
    
            #Center object around its position according to rebound
            body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
            body.sg.set_positions(center_ecl=body.center_ecl)
            self._spanglers[name]=body.sg
            
        #Set the center of the source of light for each body
        for name,body in self.bodies.items():
            body.center_source=body.source.center_ecl
            if body==self.root:
                self.center_root=body.source.center_ecl
                
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
        Misc.elapsed_time(show=True,msg="Spangling time",verbosity=VERB_NONE)
    
    def _set_observer(self,nvec=[0,0,1],alpha=0,center=None):
        """Set the position of the observer
        """
        #Only set observer if it is spangled
        if self._is_spangled():
            
            #At changing the observer, reset state
            self.sg.reset_state()
            
            #Set observer
            self.sg.set_observer(nvec=nvec,alpha=alpha,center=center)
            
            #Update areas of the spangles
            
            #Update system properties
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
            
    def _set_luz_recursive(self,name,nluz,verbosity=VERB_SIMPLE):
        """Set light source for body and 
        """
        body=self.bodies[name]
        verbose(verbosity,f"Illuminating body {name}, with nluz = {nluz} and center = {body.center_source}")
        self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
        self.sg._update_azimbeta(n_luz=nluz,name=name)
        if body.childs:
            verbose(verbosity,f"\tObject {name} has childs!")
            for child_name in body.childs:
                verbose(verbosity,f"\tCalling recursively set_luz for {child_name}")
                self._set_luz_recursive(child_name,nluz)
        else:
            verbose(verbosity,f"\tObject {name} has no childs!")
                
    def _set_luz(self,verbosity=VERB_SIMPLE):
        """Set illumination in the system.
        
        Update:
            States: illuminated, shadow, hidden_by_luz
        """
        if self._is_spangled():
            
            if not self._observer_set:
                raise AssertionError("You must first set observer before setting light.")
            
            self.bodies_illuminated=[]
            for name,body in self.bodies.items():
              
                if body.kind == "Star":
                    verbose(verbosity,f"Body {body.name} is a star... skipping")
                    continue
                    
                if body.parent.name in self.bodies_illuminated:
                    verbose(verbosity,f"Parent body of {name}, {body.parent.name}, has been already illuminated")
                    continue
                            
                #Get center of body
                center=body.center_ecl
                        
                #Get source and center
                verbose(verbosity,f"Calculating illumination for '{name}' coming from '{body.source.name}' @ {body.center_source}")            
                nluz=body.center_source-center
                        
                if body.kind == "Ring" and body.parent.kind == "Star":
                    verbose(verbosity,f"Parent body of ring, {body.parent.name} is a star. All spangles will be illuminated")
                    self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
                    cond=(self.sg.data.name==name)
                    self.sg.data.loc[cond,"unset"]=False
                    self.sg.data.loc[cond,"illuminated"]=True
                    self.sg.data.loc[cond,"shadow"]=False
                else:                
                    verbose(verbosity,f"Illuminating body {name} and all its childs")
                    self._set_luz_recursive(name,nluz,verbosity)
                    self.sg.update_illumination_state(included=body.source.shined)
                    self.bodies_illuminated+=[name]
                    
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
        """We can save time if before updating observer we save the
        illumination state

        """
        self._set_observer(nvec=self.n_obs,alpha=self.alpha_obs,center=center_obs)
        self._set_luz()
    
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
    
    def reset(self):
        """Reset system to spangling state
        """
        if self._resetable:
            self.load_from(self._snap_file_name)
            pass
        else:
            print("System is not resetable. Use resetable = True when defining the System or when you spangle it.")
    
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
        
        if self._spangled:
            
            #Integrate
            self.sim.integrate(*args,**kwargs)
            self.sim.move_to_com()
        
            #Update positions
            for name,body in self.bodies.items():
                
                #Position of the body according
                body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
    
                verbose(VERB_VERIFY,f"Updating center of body {name} @ {body.center_ecl}")
                cond=self.sg.data.name==name
                self.sg.data.loc[cond,"center_ecl"]=pd.Series([list(body.center_ecl)]*sum(cond),dtype=object).values
    
            #Update positions
            self.sg.set_positions()
            
        else:
            raise AssertionError("You must first spangle system before setting positions.")
    
    def ensamble_system(self,lamb=0,beta=0,**physics):
        """Ensamble Ringed Planet
        
        This class is for legacy purposes.
        """
        #Check if observer was provided
        if "Observer" in self.bodies:
            lamb=self.bodies["Observer"].lamb
            beta=self.bodies["Observer"].beta
            
        physics_defaults=deepcopy(LEGACY_PHYSICAL_PROPERTIES)
        physics_defaults.update(dict(limb_cs=self.bodies["Star"].limb_coeffs))
        physics_defaults.update(physics)
    
        #--CONSISTENCY--
        self._ringedplanet=dict(
            
            #Behavior
            behavior=dict(shadows=True),
            
            #Units
            CU=CanonicalUnits(UL=self.ul,UM=self.um),
    
            #Basic
            Rstar=self.bodies["Star"].radius,
            Rplanet=self.bodies["Planet"].radius,
    
            Rint=self.bodies["Ring"].fi,
            Rext=self.bodies["Ring"].fe,
            i=self.bodies["Ring"].i,
    
            a=self.bodies["Planet"].a,e=self.bodies["Planet"].e,
    
            #Orbit 
            Mstar=1,x=0,lambq=0,t0=0,kepler=False,
    
            #Observer
            eobs_ecl=np.array([lamb,beta]),
    
            #Sampling
            Np=self.bodies["Planet"].nspangles,
            Nr=self.bodies["Ring"].nspangles,
    
            Nb=0,Ns=30,
    
            #Physical properties
            physics=physics_defaults,
        )
        self.RP=RingedPlanet(**self._ringedplanet)
        return self.RP
    
