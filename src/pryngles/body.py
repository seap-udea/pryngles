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

# # Pryngles module: body 

from pryngles import *

# ## External modules

import spiceypy as spy
import numpy as np
from copy import deepcopy

# ## Aliases

sci=Science
print_df=Misc.print_df

# ## The body class
# 
# The Body class is one of the most important classes in the package. 

Body_doc="""A general body.  This class is not intended to be used independently, just for inheritance purposes.
    
Initialization attributes:

    kind : string:
        One of the kind of bodies defined in the package (see _BODY_KINDS)
        Defined objects are: "Star", "Planet", "Ring".

    defaults : OrderedDict:
        Dictionary with the properties of the object.

    primary: Class Body:
        Object in the center of the orbit of this body.

    **properties: dicitionary:
        Specification of the body properties.  All objects of the class Body has the following
        properties by default:
        
        bhash: string, default = None:
            Hash of the object, ie. a unique string identifying the object 
            (see hash Python function)

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

                prot: float [ut], default = 1:
                    Period of rotation of the star.

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

        wrot: float [rad/ut]:
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

    Create a body with None parent and hash = 'B':
    
        B=Body("Body",BODY_DEFAULTS,None,hash='B',m=2,c=2)
        
    Create a body having parent the Body "B" defined before:
         
        C=Body("Body",BODY_DEFAULTS,B,hash="C")
"""

"""
These are the default attributes for any body.
"""
BODY_DEFAULTS=dict()
BODY_DEFAULTS.update(odict(
    
    bhash=None,
    hash_by_kind=False,
    
    #Orbit
    m=1,
    x=0,y=0,z=0,
    vx=0,vy=0,vz=0,

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
))

BODY_KINDS=[]
class Body(PrynglesCommon):
    
    def __init__(self,kind,defaults,primary,**props):

        #Kind, primary and child attributes
        self.kind=kind
        self.__defaults=defaults

        #Hash object
        if 'bhash' in props:
            bhash=self.bhash=str(props["bhash"])
        elif 'hash_by_kind' in props:
            bhash=self.bhash=self.kind
        else:
            bhash=self.bhash=str(hash(self))
            

        #Update childs and parent
        if primary is not None:
            if not isinstance(primary,Body):
                raise AssertionError(f"Primary is not a valid Object: {type(primary)}, {isinstance(primary,Body)}")
            else:
                primary._update_childs(self)

        #Update primary and childs        
        self._update_primary(primary)
        self._update_childs()

        #Update default properties
        self.__dict__.update(defaults)
        #Recover hash
        self.bhash=bhash
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
                
        verbose(VERB_VERIFY,"Updating Body")
        self._update_properties()
    
    def _update_childs(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=dict()
        if child is not None:
            verbose(VERB_VERIFY,f"Add child {child.bhash} to body {self.kind} ({self.bhash})")
            self.childs[child.bhash]=child
            
    def _update_primary(self,primary=None):
        if 'primary' not in self.__dict__:
            if primary:
                verbose(VERB_VERIFY,f"Add primary {primary.bhash} to body {self.kind} ({self.bhash})")
            self.primary=primary
        elif primary is not None:
            verbose(VERB_VERIFY,f"Add parent {primary.bhash} to body {self.kind} ({self.bhash})")
            self.primary=primary
            parent._update_childs(self)
    
    def _update_properties(self):
        verbose(VERB_VERIFY,"Updating properties of Body")
        #Rotational angular velocity
        self.wrot=2*np.pi/self.prot
        #Rotation axis
        self.n_equ=sci.cartesian([1,self.roll,90*Consts.deg-self.i])

Body.__doc__=Body_doc

# ## Testing


def spangle_body(self):
    """
    Spangle the surface of the body
    """
    
    #Create spangler
    self.sg=Spangler(
        nspangles=self.nspangles,
        sphash=self.bhash,
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

    self.sg.set_observer()
    self.sg.set_luz()

Body.spangle_body=spangle_body


# ## Star Class

"""
These are the default attributes for bodies of the kind 'Star'.
"""
STAR_DEFAULTS=deepcopy(BODY_DEFAULTS)
STAR_DEFAULTS.update(odict(

    #Orbit: update
    #Same as Body
    
    #Physics: update
    #Same as Body
    
    #Optical properties: update
    limb_coeffs=[],
    spangle_type=SPANGLE_STELLAR,
    shape="sphere",
))

BODY_KINDS+=["Star"]
class Star(Body):
    """A star.

    Initialization attributes:
        
        primary: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.

            If None the object is the center of the orbit specification for other objects.
            
            Object primary for a star should be another star.
        
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
                 primary=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Star",STAR_DEFAULTS,primary,**props)

        #Check primary
        if self.primary is not None:
            if self.primary.kind!="Star":
                raise ValueError(f"Only another Star can be the primary of a Star (you provided {self.primary.kind})")

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
        self.norm_limb_darkening=Science.LIMB_NORMALIZATIONS[hash(tuple(self.limb_coeffs))]
        
    def update_star(self,**props):
        """General update propeties of the Star
        """
        verbose(VERB_VERIFY,"Updating star")
        
        Body.update_body(self,**props)
        self._update_star_properties()


# ## Planet class

"""
These are the default attributes for bodies of the kind 'Planet'.
"""
PLANET_DEFAULTS=deepcopy(BODY_DEFAULTS)
PLANET_DEFAULTS.update(odict(
    
    #Orbit: update
    #Same as Body
    
    #Physics: update
    #Same as Body
    
    #Optical: update
    spangle_type=SPANGLE_SOLID_ROCK,
    geometry="sphere",
    
    albedo_gray_spherical=1,
))

BODY_KINDS+=["Planet"]
class Planet(Body):
    """A planet.

    Initialization attributes:
        
        primary: Class Body, default = None:
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
                 primary=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Planet",PLANET_DEFAULTS,primary,**props)
        
        #Check primary
        if self.primary is None:
            raise ValueError(f"Primary not provided and it is mandatory for {self.kind}.")
        
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


# ## Ring class

RING_DEFAULTS=deepcopy(BODY_DEFAULTS)
RING_DEFAULTS.update(odict(

    #Orbit: update
    #Same as Body altough ring has not orbit properties
    
    #Physics: update
    #Same as Body
    fi=1.0,
    fe=2.0,
    
    #Optics: update
    spangle_type=SPANGLE_GRANULAR,
    shape="ring",
    albedo_gray_normal=1,
    tau_gray_optical=0,
))

# ## Ring Class

BODY_KINDS+=["Ring"]

class Ring(Body):
    """Class Ring.
    
Initialization attributes:
        
        primary: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.
            If None the object is the center of the orbit specification for other objects.

        **props: dictionary:
            List of properties for star.  For the complete set of default values of the properties
            see STAR_DEFAULTS.  Description of properties are available in the Body class documentation.
            
            Additional properties:

            fi: float [adimensional], default = 1:
                Fraction of the radius of the primary object where ring stars.

            fe: float [adimensional], default = 1:
                Fraction of the radius of the primary object where ring ends.

            albedo_gray_normal: float. default = 1: 
                Lambertian (normal) gray (wavelength indpendent) albedo of the spangle.

            tau_gray_optical: float. default = 0:
                Gray (wavelength indpendent) Optical depth of the spangle.  
                If 0 the spangle is entirely opaque to all wavelength, despite its type.            

    Derived attributes:
    
        ri: float:
            Radius of the inner border of the ring in units of the primary radius.

        re: float:
            Radius of the outer border of the ring in units of the primary radius.
            
    Notes:

        See Body class documentation.
    """
    def __init__(self,
                 primary=None,
                 **props
                ):
        
        
        #Instantiate object with basic properties
        Body.__init__(self,"Ring",RING_DEFAULTS,primary,**props)
        
        #Check primary
        if self.primary is None:
            raise ValueError(f"Primary not provided and mandatory for {self.kind}.")
        
        #Update properties
        self.update_ring(**props)

    def _update_ring_properties(self):
        """Update specific properties of the star
        
        Properties to update:
        
            ri, re: float:
                Radius of the inner (outer) border of the ring in units of the primary radius.
                Requires: limb_coefs.
                
            radius: float:
                Object radius.
                
            geometry_args: dictionary:
                
        """
        verbose(VERB_VERIFY,"Updating Ring properties")
    
        #Update radius
        self.ri=self.fi*self.primary.radius
        self.re=self.fe*self.primary.radius
        self.radius=self.re
        
        #Update geometry args for spangling purposes
        self.geometry_args=dict(ri=self.ri/self.re)
        
    def update_ring(self,**pars):
        verbose(VERB_VERIFY,"Updating Ring")
        Body.update_body(self,**pars)
        self._update_ring_properties()   


