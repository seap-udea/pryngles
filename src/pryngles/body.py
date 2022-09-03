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

#Aliases
sci=Science
print_df=Misc.print_df

# ## External modules

import spiceypy as spy
import numpy as np

# ## The body class
# 
# The Body class is one of the most important classes in the package. 

Body_doc="""A general body.  This calss is not intended to be used independently, just for inheritance purposes.
    
Initialization attributes:

    kind : string:
        One of the kind of bodies defined in the package (see _BODY_KINDS)
        Defined objects are: "Star", "Planet", "Ring".

    defaults : OrderedDict:
        Dictionary with the properties of the object.

    primary: Body
        Object in the center of the orbit of this body.

    **properties: dicitionary:
        Specification of the body properties.  All objects of the class Body has the following
        properties by default:
        
        hash: string, default = None:
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
                    Inclination of the ring with respect to the ecliptic plane.

                roll: float [rad], default = 0:
                    Roll angle.  This is the angle with respect to ecliptic x-axis in which 
                    the normal to the ring plane is rotated.

                alpha_equ: float [rad], default = 0:
                    Longitude of the zero meridian of the object.

                t0: float [ut], default = 0:
                    Initial time for zero meridian.

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
    
    hash=None,
    
    #Orbit
    m=1,

    #Physics
    radius=1,
    prot=1,
    i=0, #Inclination of the rotational axis
    roll=0,
    alpha=0, #Zero meridian
    t0=0,
    
    #Optics
    nspangles=1000,
    spangle_type=SOLID_SPANGLE,
    geometry="sphere",
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
        if 'hash' in props:
            bhash=self.hash=str(props["hash"])
        else:
            bhash=self.hash=str(hash(self))

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
        self.hash=bhash
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
                print(f"Property {prop} not identified in object {self.kind}")
                
        verbose(VERB_VERIFY,"Updating Body")
        self._update_properties()
    
    def _update_childs(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=dict()
        if child is not None:
            verbose(VERB_VERIFY,f"Add child {child.hash} to body {self.kind} ({self.hash})")
            self.childs[child.hash]=child
            
    def _update_primary(self,primary=None):
        if 'primary' not in self.__dict__:
            if primary:
                verbose(VERB_VERIFY,f"Add primary {primary.hash} to body {self.kind} ({self.hash})")
            self.primary=primary
        elif primary is not None:
            verbose(VERB_VERIFY,f"Add parent {primary.hash} to body {self.kind} ({self.hash})")
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
    self.sp=Spangler(
        nspangles=self.nspangles,
        body_hash=self.hash,
        spangle_type=self.spangle_type,
        n_equ=self.n_equ,
        alpha_equ=self.alpha,
        w_equ=self.wrot,
        t0_equ=self.t0,
    )
    
    #Populate spangler
    self.sp.populate_spangler(
        scale=self.radius,
        seed=self.seed,
        geometry=self.geometry,
        spangle_type=self.spangle_type,
        preset=self.preset,
        **self.geometry_args,
    )

Body.spangle_body=spangle_body


# ## Star Class

"""
These are the default attributes for bodies of the kind 'Star'.
"""
STAR_DEFAULTS=deepcopy(BODY_DEFAULTS)
STAR_DEFAULTS.update(odict(

    #Orbit
    
    #Physics
    
    #Optical properties
    limb_coeffs=[],
    spangle_type=STELLAR_SPANGLE,
    geometry="sphere",
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
                raise ValueError(f"Only another Star can be the primary of a Star")

        self._update_star_properties()
        
    def _update_star_properties(self):
        verbose(VERB_VERIFY,"Updating properties of Star")

        #Update limb darkening normalization
        sci.limb_darkening(0,self.limb_coeffs)
        self.norm_limb_darkening=LIMB_NORMALIZATIONS[hash(tuple(self.limb_coeffs))]
        
    def update_star(self,**props):
        verbose(VERB_VERIFY,"Updating star")
        
        Body.update_body(self,**props)
        self._update_star_properties()


# ## Planet class

"""
These are the default attributes for bodies of the kind 'Planet'.
"""
PLANET_DEFAULTS=deepcopy(BODY_DEFAULTS)
PLANET_DEFAULTS.update(odict(
    
    #Orbit
    
    #Physics
    
    #Optical
    spangle_type=SOLID_SPANGLE,
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
            
                a: float [ul], default = 1.0
                    Semi major axis of the orbit with respect to primary.

                e: float, default = 0.0
                    Eccentricity of the orbit with respect to primary.
        
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
        verbose(VERB_VERIFY,"Updating Planet properties")

    def update_planet(self,**pars):
        verbose(VERB_VERIFY,"Updating Planet")
        Body.update_body(self,**pars)
        self._update_planet_properties()


# ## Ring class

RING_DEFAULTS=deepcopy(BODY_DEFAULTS)
RING_DEFAULTS.update(odict(
    #Physics
    fi=1.0,
    fe=2.0,
    
    #Optics
    spangle_type=GRANULAR_SPANGLE,
    geometry="ring",
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
    
        ri: float [rlu]:
            Radius of the inner border of the ring

        re: float [rlu]:
            Radius of the outer border of the ring
            
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


