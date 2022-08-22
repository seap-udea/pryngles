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

# # Pryngles module: stars

from pryngles import *

# ## External modules

sci=Science
print_df=Misc.print_df

# ## Star Class

"""
These are the default attributes for bodies of the kind 'Star'.

DEVELOPER:
    You may add new attributes as the model gets more complex.
    Please document properly each attribute.

"""
STAR_DEFAULTS=odict(

    #Orbit
    m=1,
    
    #Physics
    radius=1,
    prot=1,
    i=0, #Inclination of the rotational axis
    roll=0,
    alpha=0, #Zero meridian
    t0=0,
    
    #Optical properties
    nspangles=1000,
    limb_coeffs=[],
)

BODY_KINDS+=["Star"]
class Star(Body):
    """A star.

    Initialization attributes:
        
        primary: Class Body, default = None:
            Object in the center of the orbit of the star for specification purposes.

            If None the object is the center of the orbit specification for other objects.
            
            Object primary for a star should be another star.
        
        orbit:

            These attributes should be compatible with rebound.

            m: float [um], default = 1: 
                Mass of the star. It should be different than zero.

        physics:

            radius: float [ul], default = 1:
                Radius of the star.

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

        optics:

            nspangles: int, default = 1000:
                Number of spangles on which the star will be discretized.

            limb_coeffs: list [adimensional], default = []
                List of limb darkening fit coefficients.  See Science.calc_limbdarkening.
                
                Models in: https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html
                Coefficients available at: https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt

    Derived attributes:
    
            wrot: float [rad/ut]:
                Rotational angular velocity.
                
            n_equ: array(3):
                Rotational axis vector.
    
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
        Body.__init__(self,"Star",STAR_DEFAULTS,primary)

        #Check primary
        if self.primary is not None:
            if self.primary.kind!="Star":
                raise ValueError(f"Only another Star can be the primary of a Star")
                
        #Update properties
        self.update_body(**props)
        
    def update_body(self,**props):
        Body.update_body(self,**props)
        
        #Update physics
        
        #Rotational angular velocity
        self.wrot=2*np.pi/self.prot
        
        #Rotation axis
        self.n_equ=sci.cartesian([1,self.roll,90*Consts.deg-self.i])


def spangle_body(self,seed=0,preset=False):
    """
    Spangle the surface of the star
    """
    
    #Create spangler
    self.sp=Spangler(
        nspangles=self.optics.nspangles,
        body_hash=self.hash,
        spangle_type=STAR_SPANGLE,
        n_equ=self.physics.n_equ,
        alpha_equ=self.physics.alpha,
        w_equ=self.physics.wrot,
        t0_equ=self.physics.t0,
    )
    
    #Populate spangler
    self.sp.populate_spangler(
        scale=self.physics.radius,
        seed=seed,
        geometry="sphere",
        preset=preset
    )

Star.spangle_body=spangle_body


