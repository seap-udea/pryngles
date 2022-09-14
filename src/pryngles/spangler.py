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

# # Pryngles module: Spangler

# This module contains all the physics of light scattered on spangles

from pryngles import *

# ## External modules

import pandas as pd
import random

#Specialized plotting methods
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib import animation
from celluloid import Camera # getting the camera
from ipywidgets import interact,fixed,widgets

# ## Aliases

print_df=Misc.print_df
sci=Science

# ## Constants

#Colors: Given in hue (0-360), level (0: black-1: white), saturation (0-1)

#Type of spangles
SPANGLE_COLORS=dict()

#Planetary spangles
SOLID_SPANGLE=0
SPANGLE_COLORS[SOLID_SPANGLE]=[27,0.5,1.0]

ATMOSPHERIC_SPANGLE=1
SPANGLE_COLORS[ATMOSPHERIC_SPANGLE]=[27,0.5,1.0]

LIQUID_SPANGLE=2
SPANGLE_COLORS[LIQUID_SPANGLE]=[195,0.7,0.5]

#Ring or disks spangles
GRANULAR_SPANGLE=3
SPANGLE_COLORS[GRANULAR_SPANGLE]=[59,0.7,0.4]

#Semitransparent
GASEOUS_SPANGLE=4
SPANGLE_COLORS[GASEOUS_SPANGLE]=[27,0.5,1.0]

#Stellar spangles
STELLAR_SPANGLE=5
SPANGLE_COLORS[STELLAR_SPANGLE]=[59,0.7,1.0]

#List of semitransparent spangles
SEMITRANSPARENT_SPANGLES=[GRANULAR_SPANGLE,GASEOUS_SPANGLE]

#Color of shadow
SHADOW_COLOR=[180,0.6,0.0]
SHADOW_COLOR_OBS=[180,0.2,0.0]

# ## The Spangler class
# 
# This class contains a family of routines useful for spangling different kind of objects.

Spangler_doc="""A Spangler associated to an object or set of objects.
    
   There are two ways to initialize a Spangler:
    
        Creating a Spangler for a single object:
        
            Mandatory:

                nspangles: int, default = 0:
                    Number of spangles in spangling.

            Optional:

                body_hash: string, default = None:
                    Hash identifying the body to which spangles are associated 
                    (see Body documentation for explanation about hash).

                spangle_type: int, default = 0:
                    Type of spangle (see *_SPANGLE in Consts module).

                n_equ: numpy Array (3), default = [0,0,1]:
                    unitary vector normal to {equ} (equatorial) plane.

                alpha_equ: float, default = 0:
                    Roll angle of x-axis of equatorial system (not implemented yet)

                center_equ: numpy Array (3), default = [0,0,0]:
                    Position of the spnagler in the {equ} (equatorial) system.

                center_ecl: numpy Array (3), default = [0,0,0]:
                    Position of the spnagler in the {ecl} (ecliptic) system.
                    
                w, q0: float [rad/ut, rad], default = 0, 0:
                    Angular velocity and reference latitude at t = 0.

        Joining a set of Spanglers (several objects):

            spanglers: list of Spanglers. default = []:
                Set of spanglers to join.

Core attributes:

    nspangles: int:
        Total number of spangles.

    data: Pandas DataFrame: 
        Dataframe containing all the information about the spangling.
        For Columns see global variable SPANGLER_COLUMNS.
"""

#Columns of spangling
SPANGLER_COLUMNS=odict({
    "sphash":"",

    #Type of spangle
    "spangle_type":SOLID_SPANGLE, #For a list of spangle types see the constants module.
    "dim":2, #Dimension of the spangle, 3: sphere, 2: circle, ring

    #Lengh-scale
    "scale":1, #The length scale of the body, eg. for a ring this is the outer radius

    #Body parameters
    "n_equ":[0,0,1], #Direction of the equator of the body with respect
    "alpha_equ":0, #Zero meridian of equatorial system
    "w":0, #Rotational angular velocity [rad/ut]
    "q0":0, #Initial time [rad], Longitude (azimutal angle) are calculated as: q = q0 + w (t - t0)

    #Coordinates of the spangle (cartesian and spherical) in the body-centric system
    "center_equ":[0,0,0],#Center of the body with respect to barycenter
    "x_equ":1,"y_equ":0,"z_equ":0, #Cartesian coordinates
    "r_equ":1,"q_equ":0,"f_equ":0, #Spherical coordinates: q: longitude, f: latitude
    "ns_equ":[0,0,1], #Unitary vector normal to the spangle

    #Coordinates of the spangle (cartesian and spherical) in the ecliptic system
    "center_ecl":[0,0,0], #Center of the body with respect to barycenter
    "x_ecl":1,"y_ecl":0,"z_ecl":0, #Cartesian coordinates of the spangle
    "ns_ecl":[0,0,1],#Unitary vector normal to the spangle, calculated in the class

    #Coordinates of the spangle (cartesian and spherical) in the observer system
    "center_obs":[0,0,0], #Center of the body with respect to observer
    "x_obs":1,"y_obs":0,"z_obs":0, #Cartesian coordinates of the spangle
    "ns_obs":[0,0,1],#Unitary vector normal to the spangle, calculated in the class
    "rho_obs":1,"az_obs":0,"cost_obs":0, #Cylindrical coordinates of the spangle: rho, q, cos(theta)
    "cos_obs":1, #Angle between normal to spangle and direction of observer
    "d_obs":1, #Distance of the Spangle to light-source
    
    #Coordinates of the spangle (cartesian and spherical) in the light-source system
    "center_luz":[0,0,0],#Center of the body with respect to light
    "x_luz":1,"y_luz":0,"z_luz":0,#Calculated in the class
    "ns_luz":[0,0,1],#Unitary vector normal to the spangle, calculated in the class
    "rho_luz":1,"az_luz":0,"cost_luz":0, #Cylindrical coordinates of the spangle: rho, q, cos(theta)
    "cos_luz":1, #Angle between normal to spangle and direction of light-source
    "d_luz":1, #Distance of the Spangle to light-source

    #Coordinates of the spangle (cartesian and spherical) in the intersection system
    "center_int":[0,0,0],#Center of the body 
    "x_int":1,"y_int":0,"z_int":0,#Cartesian coordinates
    "ns_int":[0,0,1],#Unitary vector normal to the spangle, calculated in the class
    "rho_int":1,"az_int":0,"cost_int":0, #Pseudo cylindrical coordinates (rho, q, cos(theta))
    "cos_int":1, #Angle between normal to spangle and direction of intersection
    "d_int":1, #Distance of the Spangle to intersection

    #Geometrical parameters
    "asp":1.0, #Effective area of the spangle
    "dsp":1.0, #Effective diameter of spangle, dsp = 2*(asp/pi)**0.5

    #Optical parameters
    "albedo_gray_normal":1.0,
    "tau_gray_optical":0.0,
    
    #Special states
    "unset":True, #State has not been set
    "hidden":False, #The spangle is not taken into account for photometry
    "source":False, #The spangle belongs to a light-source (it does not reflect light)
})

SPANGLER_VISIBILITY_STATES=odict({
    #Spangle state
    "visible":False, #The spangle is visible from observer
    "shadow":False, #The spangle is in the shadow of other spangle
    "indirect":False, #The spangle is indirectly illuminated
    "emit":False, #The spangle is emmitting
    "intersect":False, #Intermediate state to calculate intersections
    "above":False, #Intermediate state to calculate above or below state respect to ring
})
SPANGLER_COLUMNS.update(SPANGLER_VISIBILITY_STATES)

SPANGLER_SOURCE_STATES=odict({
    "illuminated":False, #The spangle is illuminated by the light-source
    "transmit":False, #The spangle is illuminated but transmitting light
    "transit":False, #The spangle is transiting
    "occult":False, #The spangle is occulted by a light source
})
SPANGLER_COLUMNS.update(SPANGLER_SOURCE_STATES)

SPANGLER_COL_COPY=["center","x","y","z","ns","rho","az","cost","cos","d"]
SPANGLER_COL_LUZ=[column+"_luz" for column in SPANGLER_COL_COPY]
SPANGLER_COL_OBS=[column+"_obs" for column in SPANGLER_COL_COPY]
SPANGLER_COL_INT=[column+"_int" for column in SPANGLER_COL_COPY]

SPANGLER_RING_BORDER=100

class Spangler(PrynglesCommon):
    
    def __init__(self,
                 #Initialization using specific options
                 #Basic
                     nspangles=1,
                     sphash=None,
                     n_equ=SPANGLER_COLUMNS["n_equ"],
                     alpha_equ=SPANGLER_COLUMNS["alpha_equ"],
                     center_equ=SPANGLER_COLUMNS["center_equ"],
                     center_ecl=SPANGLER_COLUMNS["center_ecl"],
                 #Optional
                     w=SPANGLER_COLUMNS["w"],
                     q0=SPANGLER_COLUMNS["q0"],
                 #Initialization with a list of 
                     spanglers=[]
                ):
        
        #Common attributes
        self.n_obs=[0,0,1]
        self.n_luz=[0,0,1]
        #Direction of the observer in spherical coordinates
        self.rqf_obs=[0,0,90*Consts.deg]
        #Transformation matrices from equatorial to ecliptic coordinates
        self.M_equ2ecl=dict()
        self.M_ecl2equ=dict()
        #Convex hulls of spanglers
        self.qhulls=dict()
        
        #Create a spanglers with a list of other spanglers
        if len(spanglers)>0:
            verbose(VERB_SIMPLE,f"Joining {len(spanglers)} spanglers")
            self._join_spanglers(spanglers)
            
        #Create a spangler with the desired options
        else:
            #Attributes
            self.nspangles=nspangles
            self.geometry="Vanilla"
            
            #Default property values
            self._defaults=deepcopy(SPANGLER_COLUMNS)

            if not sphash:
                self.sphash=str(random.getrandbits(16))
                verbose(VERB_VERIFY,f"Generating random hash {self.sphash}")
            else:
                self.sphash=sphash
                
            self._defaults.update(dict(sphash=self.sphash))
            
            #Update other parameters
            self._defaults.update(
                dict(w=w,q0=q0)
            )

            #Create Spangler dataframe
            if self.nspangles>0:
                
                #Create a simple DataFrame with the default values
                self.data=pd.DataFrame([list(self._defaults.values())]*self.nspangles,
                                       columns=self._defaults.keys())

                #Update positions
                self.set_positions(
                    n_equ=n_equ,alpha_equ=alpha_equ,
                    center_equ=center_equ,center_ecl=center_ecl,
                    t=None
                )
        
            else:        
                verbose(VERB_SIMPLE,f"Creating a blank Spangler")
                #Creat a blank DataFrame
                self.data=pd.DataFrame(columns=self._defaults.keys())
        
    def set_rotation(self,sphash,w,t0):
        """Set rotational parameters
        """
        cond=(self.data.sphash==sphash)
        self.data.loc[cond,"w"]=w
        self.data.loc[cond,"t0"]=t0
    
    def reset_state(self):
        """Reset spangler state
        """
        self.data[list(SPANGLER_SOURCE_STATES)+list(SPANGLER_VISIBILITY_STATES)]=False
        self.data["unset"]=True

    def set_scale(self,scale):
        """Set scale
        """
        lengths=[
            "x_equ","y_equ","z_equ",
            "x_ecl","y_ecl","z_ecl",
            "x_obs","y_obs","z_obs","d_obs",
            "x_luz","y_luz","z_luz","d_luz",
            "r_equ","rho_obs","rho_luz",
            "dsp",
        ]
        self.scale=scale
        self.data[lengths]*=self.scale
        areas=[
            "asp",
        ]
        self.data[areas]*=self.scale**2
        vectors=[
            "center_ecl",
            "center_equ",
            "center_obs",
        ]
        for vector in vectors:
            self.data[vector]=[np.array(v)*scale for v in self.data[vector]]
        
    def _join_spanglers(self,spanglers):
        """
        Join spanglers into a single spangler

        Parameters:
            spanglers: list of Spanglers:
                Spanglers to join.
        """
        self.sphash=[]
        for spangler in spanglers:
            if not isinstance(spangler,Spangler):
                raise AssertionError(f"One of the spangler is not an Spangler instance")
                
            if spangler.sphash in self.sphash:
                raise ValueError(f"Hash '{spangler.sphash}' already included in spangler '{self.sphash}'")
                
            self.sphash+=[spangler.sphash]
        
        #Set of spanglers
        self.spanglers=spanglers

        #Concatenate data
        datas=[spangler.data for spangler in spanglers]
        self.data=pd.concat(datas,ignore_index=True)

        self.M_equ2ecl=dict()
        for spangler in spanglers:
            self.M_equ2ecl.update(spangler.M_equ2ecl)

        #Join properties
        self.nspangles=len(self.data)
        self.geometry="Join"
    
Spangler.__doc__=Spangler_doc

# ### Set positions

def set_positions(self,
                  n_equ=[],alpha_equ=0,
                  center_equ=[],center_ecl=[],
                  t=None
                 ):
    """
    Set the positions and orientation of spanglers in all reference systems.

    Parameters:

        n_equ: list/array (3), default = []:
            Normal vector towards north pole equatorial system.

        alpha_equ: float, default = 0:
            Roll angle of x-axis of equatorial system (not implemented yet)

        center_equ: list/array (3), default = []:
            Location of the center of the body with respect to the barycenter.
            
        center_ecl: list/array (3), default = []:
            Location of the center of the body with respect to the barycenter.
            
        t: float, default = None:
            Time.  This quantity is used to update the equatorial coordinates.
            If None, equatorial coordinates are not set.

    Return:
        None

    Update:

        Coordinates of the spangles, (x_ecl,y_ecl,z_ecl).

        Rotation matrices M_equ2ecl
    """
    verbose(VERB_VERIFY,f"Setting positions")

    #Update normal vectors
    qupdate=False
    
    #Update center
    if len(center_equ)>0:
        verbose(VERB_VERIFY,f"Updating center in {{equ}} to {center_equ}")
        self.data["center_equ"]=[center_equ]*self.nspangles
    if len(center_ecl)>0:
        verbose(VERB_VERIFY,f"Updating center {{ecl}} to {center_ecl}")
        self.data["center_ecl"]=[center_ecl]*self.nspangles
        
    if len(n_equ)>0:
        verbose(VERB_VERIFY,f"Generating equatorial transformation matrices from n_equ = {n_equ}")
        
        #Unitary equatorial vector
        n_equ,one=spy.unorm(n_equ)
        self.data["n_equ"]=[n_equ]*self.nspangles

        #Transformation matrices
        self.M_equ2ecl[self.sphash],M_ecl2equ=sci.rotation_matrix(n_equ,alpha_equ)
        
        qupdate=True

    #Update equatorial coordinates by rotation
    if t is not None:
        verbose(VERB_VERIFY,f"Updating rotations at t = {t}")

        self.data["q_equ"]=[q+q0+w*t for q,w,q0 in zip(self.data.q_equ,self.data.w,self.data.q0)]
        self.data[["x_equ","y_equ","z_equ"]]=            [sci.cartesian(r) for r in np.array(self.data[["r_equ","q_equ","f_equ"]])]

        qupdate=True
        
    if qupdate:
        #Update normal vectors
        cond=(self.data.dim==2)
        self.data.loc[cond,"ns_equ"]=pd.Series([[0,0,1]]*cond.sum(),dtype=object)
        lista=[spy.unorm(list(r))[0] for r in np.array(self.data[~cond][["x_equ","y_equ","z_equ"]])]
        self.data.loc[~cond,"ns_equ"]=pd.Series(lista,dtype=object)
        
    #Convert from equatorial to ecliptic
    verbose(VERB_VERIFY,f"Converting to equatorial")
    self.data[["x_ecl","y_ecl","z_ecl"]]=        [np.matmul(self.M_equ2ecl[sph],r+cequ)+cecl         for sph,r,cequ,cecl in zip(self.data.sphash,
                                   np.array(self.data[["x_equ","y_equ","z_equ"]]),
                                   self.data.center_equ,self.data.center_ecl)]
    
    #Update spangles orientations
    verbose(VERB_VERIFY,f"Generating normal vectors")
    self.data["ns_ecl"]=[np.matmul(self.M_equ2ecl[sph],n) for sph,n in zip(self.data.sphash,
                                                                         self.data.ns_equ)]
    
    #Update velocities
    #Not implemented yet
    
Spangler.set_positions=set_positions

# ### Test class




# ### Populate Spangler

def populate_spangler(self,
                      preset=False,
                      scale=1,seed=0,spangle_type=SOLID_SPANGLE,
                      geometry="circle",**geometry_args):
    
    """Populate data of a Spangler using points generated with a given geometry.
    
    Parameters:
    
        preset: boolean, default = False:
            If true the spangler is populated with preset data.
            
        scale: float. default = 1:
            Scale size of the object.
            
        seed: integer. default = 0:
            Value of the integer seed of random number generation (if 0 no random seed is set).
            If a non-zero seed is used the position of the spangle will be always the same.
            
        geometry: string, default = "circle":
            Geometry of the Sampler.  Available: "circle", "ring", "sphere"
            
        geometry_args: dictionary:
            See Sampler methods documentation.
             
    """
    #Check if preset
    if preset:
        verbose(VERB_VERIFY,f"Populating spangler from preset for {geometry}")
        preset=(geometry,geometry_args)
        self.sample=Sampler(preset=preset,N=self.nspangles,seed=seed)   
    else:
        self.sample=Sampler(N=self.nspangles,seed=seed)
        exec(f"self.sample.gen_{geometry}(**geometry_args)")

    self.geometry=geometry
    self.data["dim"]=self.sample.dim
    self.data["spangle_type"]=spangle_type

    if self.sample.dim>2:
        #Purge sample if it is in 3d
        verbose(VERB_VERIFY,f"Purging 3d sample")
        self.sample.purge_sample()
        
    elif self.geometry == "ring":
        #Add hidden spangles to ring inner borders
        pp_border=np.zeros((SPANGLER_RING_BORDER,3))
        ss_border=np.zeros((SPANGLER_RING_BORDER,3))
        for i,theta in enumerate(np.linspace(0,2*np.pi,SPANGLER_RING_BORDER)):
            pp_border[i]=[self.sample.ri,theta,0]
            ss_border[i]=[self.sample.ri*np.cos(theta),
                          self.sample.ri*np.sin(theta),
                          0]
        self.sample.pp=np.vstack((self.sample.pp,pp_border))
        self.sample.ss=np.vstack((self.sample.ss,ss_border))
        self.sample.N+=SPANGLER_RING_BORDER
                
    #Check if number of samples is not equal to that of spangles
    if self.sample.N!=self.nspangles:
        verbose(VERB_SYSTEM,f"Sample size {self.sample.N} is different from spangles {self.nspangles}. Adjusting.")
        dif=self.sample.N-self.nspangles
        if dif>0:
            verbose(VERB_SYSTEM,f"Adding {dif} entries to DataFrame")
            for i in range(dif):
                df=pd.DataFrame([self.data.iloc[-1]])
                self.data=pd.concat([self.data,df],ignore_index=True)
        else:
            verbose(VERB_SYSTEM,f"Removing {-dif} entries to DataFrame")
            self.data.drop(range(self.nspangles+dif,self.nspangles),inplace=True)
        self.nspangles=self.sample.N
    
    #Area
    self.data["asp"]=self.sample.aes*scale**2
    self.data["dsp"]=2*(self.data["asp"]/np.pi)**0.5
    
    #Update scale
    self.data["scale"]=scale

    #Store positions in DataFrame
    self.data[["x_equ","y_equ","z_equ"]]=self.sample.ss*scale
    self.data[["r_equ","q_equ","f_equ"]]=self.sample.pp
    self.data["r_equ"]*=scale

    #Update normal vectors
    if self.sample.dim>2:
        self.data["ns_equ"]=[spy.unorm(list(r))[0] for r in np.array(self.data[["x_equ","y_equ","z_equ"]])]
    else:
        self.data["ns_equ"]=pd.Series([[0,0,1]]*self.nspangles,dtype=object)
        
    #Hide border points in case of ring
    if geometry == "ring":
        self.data.loc[self.nspangles-SPANGLER_RING_BORDER:self.nspangles,"hidden"]=True
        
    #Update positions
    self.set_positions()
    
Spangler.populate_spangler=populate_spangler


# ### Plot3D

def plot3d(self,
           center_at=None,
           not_plot=[],
           coords="ecl",
           factor=1.2,
           show_hidden=True,
           figsize=(5,5)
          ):
    """Plot spangle in 3d.

    Optional parameters:
    
        center_at: string, default = None:
            Hash of the object around which the plotting will be centered at (see sphash column
            of the Spangler DataFrame).
            
        not_plot: list of strings, default = []:
            List of object hashes to not plot.
            
        coords: list of strings, default = ["x_ecl","y_ecl","z_ecl"]:
            Which coordinates do you want to plot.  
            Available: equ, ecl, obs, luz, int.

        factor: float, default = 1.2:
            Size of the coordinate axes.  factor = 1 correspond to axis equal to maximum and minumum.

        show_hidden: boolean, default = True:
            If True show hidden spangles (used to create convex hull of objects).
            
    """
    bgcolor='k'

    #Check if plot is in the ecliptic system
    qecl=True
    if 'ecl' not in coords:
        qecl=False
    scoords=coords
    coords=[f"x_{scoords}",f"y_{scoords}",f"z_{scoords}"]
    
    #Center
    cond=(self.data.sphash==center_at)
    x_cen,y_cen,z_cen=self.data[cond][coords].mean() if sum(cond)>0 else np.array([0,0,0])
    
    #Figure
    fig=plt.figure(figsize=figsize)
    fig.patch.set_facecolor(bgcolor)
    ax=fig.add_subplot(111,projection='3d',facecolor=bgcolor)
    ax.axis("off")
    
    #Spangles
    for i in range(self.nspangles):

        #Avoid plotting 
        sphash=self.data.loc[i,"sphash"]
        if sphash in not_plot:
            continue
        
        #Reference transparency of spangles
        alpha_base=0.5

        #Avoid hidden spangles
        if self.data.loc[i,"hidden"]:
            continue

        spangle_type=self.data.loc[i,"spangle_type"]

        #Define color according to illumination
        if self.data.loc[i,"illuminated"]:
            color=Misc.rgb(SPANGLE_COLORS[spangle_type]) #Planet color
        else:
            color=Misc.rgb(SHADOW_COLOR) #Gray

        if self.data.loc[i,"transmit"]:
            alpha_base=0.2

        #Define alpha according to albedo
        alpha=alpha_base*self.data.albedo_gray_normal[i]

        center=[self.data[coords[0]][i]-x_cen,self.data[coords[1]][i]-y_cen,self.data[coords[2]][i]-z_cen]
        radius=self.data.dsp[i]/2
        zDir=self.data[f"ns_{scoords}"][i]

        #Change properties according to visibility
        ec=color
        lw=0
        if ~self.data.loc[i,"visible"]:
            ec='k'
            lw=0
            alpha*=0.5

        #verbose(VERB_DEEP,i,center,radius,zDir)
        Plot.circle3d(ax,
                      center=center,
                      radius=radius,
                      zDir=zDir,
                      color=color,alpha=alpha,ec=ec,lw=lw)

    #Scatter plot of transmit
    cond=(~self.data.hidden)&(self.data.transmit)
    ax.scatter(self.data[cond][coords[0]],self.data[cond][coords[1]],self.data[cond][coords[2]],
               marker='v',s=10,ec='w',fc='None',alpha=0.1)
        
        
    #Aspect
    ax.set_box_aspect([1,1,1])

    #Zoom around center
    cond=(self.data.sphash==center_at)
    cond=cond if sum(cond)>0 else [True]*self.nspangles

    #Not 
    cond=cond&(~self.data.sphash.isin(not_plot))
    
    #Range
    maxval=1.0*np.abs(self.data[cond][coords].to_numpy()-[x_cen,y_cen,z_cen]).max()
    ax.set_xlim(-maxval,maxval)
    ax.set_ylim(-maxval,maxval)
    ax.set_zlim(-maxval,maxval)
    
    #Decoration
    xmin,xmax=factor*np.array(list(ax.get_xlim()))
    ymin,ymax=factor*np.array(list(ax.get_ylim()))
    zmin,zmax=factor*np.array(list(ax.get_zlim()))

    #Axis
    ax.plot([xmin,xmax],[0,0],[0,0],'w-',alpha=0.3)
    ax.plot([0,0],[ymin,ymax],[0,0],'w-',alpha=0.3)
    ax.plot([0,0],[0,0],[zmin,zmax],'w-',alpha=0.3)
    ax.text(xmax,0,0,rf"$x_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
    ax.text(0,ymax,0,rf"$y_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
    ax.text(0,0,zmax,rf"$z_{{{scoords}}}$",color='w',alpha=0.5,fontsize=8)
    
    #Plot n_obs and n_luz vector only in the case of ecliptic system
    if qecl:
        increase=1.05*factor*maxval
        if "n_luz" in self.__dict__:
            #Light
            ax.quiver(+self.n_luz[0]*increase,+self.n_luz[1]*increase,+self.n_luz[2]*increase,
                      -self.n_luz[0]*increase,-self.n_luz[1]*increase,-self.n_luz[2]*increase,
                      color='y',alpha=0.7)
            ax.text(self.n_luz[0]*increase,self.n_luz[1]*increase,self.n_luz[2]*increase,
                    r"$n_{luz}$",color='w',alpha=0.7,fontsize=8,ha='left',va='bottom')

        if "n_obs" in self.__dict__:
            #Observer
            ax.quiver(+self.n_obs[0]*increase,+self.n_obs[1]*increase,+self.n_obs[2]*increase,
                      -self.n_obs[0]*increase,-self.n_obs[1]*increase,-self.n_obs[2]*increase,
                      color='c',alpha=0.7)
            ax.text(self.n_obs[0]*increase,self.n_obs[1]*increase,self.n_obs[2]*increase,
                    r"$n_{obs}$",color='c',alpha=0.7,fontsize=8,ha='right',va='top')

            r_obs,t_obs,f_obs=sci.spherical(self.n_obs)
            ax.view_init(f_obs*Consts.rad,t_obs*Consts.rad)
    else:
        ax.view_init(30,60)
        
    #Title
    ax.set_title(f"Spangler {self.geometry}, N = {self.nspangles}",
                 color='w',fontsize=10)
    Plot.pryngles_mark(ax)
    
    #Scale
    ax.text2D(0,0,f"Axis scale: {maxval*factor:.2g}",
            fontsize=8,color='w',
            transform=ax.transAxes)

    fig.tight_layout()
    self.fig3d=fig
    self.ax3d=ax

Spangler.plot3d=plot3d


# ### Set intersection, observer, light-source

def set_intersect(self,
                  nvec,alpha=0,
                  center=np.array([0,0,0]),
                  sphash=None
                 ):
    """Set the positions and orientation of spanglers in an intersection direction

    Parameters:

        nvec: list/array (3), default = []:
            Normal vector towards the observer.

        alpha: float, default = 0:
            Roll angle of x-axis.
            
    Optional:
        
        center: list/array(3), default = None:
            Define the position of the vantage point in the ecliptic system.
            
        sphash: string, default = None:
            Spangler hash to which the transformation will be applied.
            
    Return:
        None

    Update:

        qhulls: dictionary with convex hulls of bodies.

        Coordinates of the spangles, (x_int,y_int,z_int).

        Normal to spangles, ns_int.

        Rotation matrices M_ecl2obs, M_int2ecl, 
    """
    verbose(VERB_SIMPLE,f"Setting intersect")
    
    #Ensure that center is array
    center=np.array(center)

    verbose(VERB_VERIFY,f"Generating intersect matrices from nvec = {nvec}")
    #Unitary observer vector
    n_int,d_int=spy.unorm(nvec)
    alpha_int=alpha
    rqf_int=sci.spherical(n_int)

    #Transformation matrices
    M_int2ecl,M_ecl2int=Science.rotation_matrix(n_int,alpha_int)
    
    #Depending on body
    cond=[True]*self.nspangles
    if sphash:
        cond=self.data.sphash==sphash

    #Update positions
    self.data.loc[cond,["x_int","y_int","z_int"]]=        [np.matmul(M_ecl2int,r-center) for r in np.array(self.data[cond][["x_ecl","y_ecl","z_ecl"]])]
    
    #Center of the object in the observer reference system
    self.data.loc[cond,"center_int"]=        pd.Series([np.matmul(M_ecl2int,r-center) for r in np.array(self.data.center_ecl)])
    
    #Pseudo-cylindrical coordinates in the observer system
    self.data.loc[cond,["rho_int","az_int","cost_int"]]=        [sci.pcylindrical(r) for r in np.array(self.data[cond][["x_int","y_int","z_int"]])-np.vstack(self.data[cond].center_int)]

    #Compute distance to light-source of each spangle
    self.data.loc[cond,"d_int"]=d_int-self.data[cond].z_int #Asuming d_int > radius of the object

    #Direction of spangle with respect to direction
    self.data.loc[cond,"cos_int"]=[np.dot(n_ecl,n_int) for n_ecl in self.data.ns_ecl[cond]]
    
    #Update spangles orientations
    lista=[np.matmul(M_ecl2int,n_ecl) for n_ecl in self.data[cond].ns_ecl]
    self.data.loc[cond,"ns_int"]=pd.Series(lista,dtype=object).values

    #Spherical
    self.rqf_obs=sci.spherical(nvec)
    
    #Convex hulls
    for sphash in Misc.flatten([self.sphash]):

        self.qhulls[sphash]=[]
        cond_obj=(self.data.sphash==sphash)
        
        if (self.data[cond_obj].hidden).sum()==0:
            #Convex hull of whole objects
            cond_hull=(cond_obj)&(~self.data[cond_obj].hidden)
            verbose(VERB_SIMPLE,"Hull points (whole object):",sum(cond_hull))
            zmin,zmax=self.data[cond_hull]["z_int"].values.min(),self.data[cond_hull]["z_int"].values.max()
            zmed=0.5*(zmin+zmax)
            self.qhulls[sphash]+=[dict(
                sphash=sphash,
                hulltype="med",
                zmin=zmin,zmax=zmax,zmed=zmed,
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]])
            )]
        else:
            #Convex hull of objects with a hole (eg. rings)
            
            #Plane of rings
            cond_hidden=(cond_obj)&(self.data[cond_obj].hidden)
            hidden=self.data[cond_hidden][["x_int","y_int","z_int"]].values
            nhidden=len(hidden)
            p1,p2,p3=hidden[0],hidden[int(nhidden/3)],hidden[2*int(nhidden/3)]
            plane=Science.Plane(p1,p2,p3)
            
            #Convex hull of hidden points (the hole)
            cond_hull=(cond_obj)&(self.data[cond_obj].hidden)
            verbose(VERB_SIMPLE,"Hull points (hidden):",sum(cond_hull))
            
            zmin,zmax=self.data[cond_hull]["z_int"].values.min(),self.data[cond_hull]["z_int"].values.max()
            zmed=0.5*(zmin+zmax)

            self.qhulls[sphash]+=[dict(
                sphash=sphash,
                hulltype="hidden",
                zmin=zmin,zmax=zmax,zmed=zmed,
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]]),
                plane=plane
            )]
        
            #Convex hull of no hidden points
            cond_hull=(cond_obj)&(~self.data[cond_obj].hidden)

            #Normal vector to ring
            ns_int=self.data[cond_hull]["ns_int"].iloc[0]

            verbose(VERB_SIMPLE,"Hull points (visible ring):",sum(cond_hull))
            
            zmin,zmax=self.data[cond_hull]["z_int"].values.min(),self.data[cond_hull]["z_int"].values.max()
            zmed=0.5*(zmin+zmax)
            self.qhulls[sphash]+=[dict(
                sphash=sphash,
                hulltype="plane",
                zmin=zmin,zmax=zmax,zmed=zmed,
                ns_int=ns_int,
                qhull=Science.get_convexhull(self.data[cond_hull][["x_int","y_int"]]),
                plane=plane
            )]
            
    return cond,n_int,d_int

Spangler.set_intersect=set_intersect


def set_observer(self,nvec=[0,0,1],alpha=0,center=np.array([0,0,0])):
    """Set the positions and orientation of spanglers in the observer system.

    Parameters:

        nvec: list/array (3), default = []:
            Normal vector towards the observer.

        alpha: float, default = 0:
            Roll angle of x-axis of observer system (not implemented yet)
    """
    verbose(VERB_SIMPLE,f"Setting observer")
    cond,self.n_obs,self.d_obs=self.set_intersect(nvec,alpha,center)
    self.data.loc[cond,SPANGLER_COL_OBS]=self.data.loc[cond,SPANGLER_COL_INT].values

Spangler.set_observer=set_observer


def set_luz(self,nvec=[0,0,1],center=np.array([0,0,0]),sphash=None):
    """Set the positions and orientation of spanglers in the light-source system.

    Parameters:

        nvec: list/array (3), default = []:
            Normal vector towards the observer.

        sphash: string, default = None:
            Body to apply this light direction

    """
    verbose(VERB_SIMPLE,f"Setting light-source")
    cond,self.n_luz,self.d_luz=self.set_intersect(nvec=nvec,center=center,sphash=sphash)
    self.data.loc[cond,SPANGLER_COL_LUZ]=self.data.loc[cond,SPANGLER_COL_INT].values

Spangler.set_luz=set_luz


# ### Simple state updater

def update_simple_state(self):
    """Update spangle states with the simplest algorithm
    """
    self.data.unset=False
    
    #Condition for visibility
    """
    & ! Hidden
        (
            | cos_obs > 0: spangler it is towards the observer
            | Spangle type is semitransparent
        )
    """
    cond=    (~self.data.hidden)&    (        (self.data.cos_obs>0)|        (self.data.spangle_type.isin(SEMITRANSPARENT_SPANGLES))
    )
    self.data.loc[cond,"visible"]=True
    
    #Condition for illumination
    """
    & ! Hidden
        (
            | dim = 2: 2d spangles are always illuminated
            | cos_luz > 0: spangler it is towards the light source
            | Spangle type is stellar
        )
    """
    cond=    (~self.data.hidden)&    (        (self.data.dim==2)|        (self.data.cos_luz>0)|        (self.data.spangle_type==STELLAR_SPANGLE)
    )
    self.data.loc[cond,"illuminated"]=True

    #Conditions for transmission:
    """
    & No hidden
    (
        & Spangle type is semitransparent
        & cos_obs . cos_luz < 0: observer and light source are in opposite sides
    )
    """
    #"""
    cond=    (~self.data.hidden)&    (     (self.data.spangle_type.isin(SEMITRANSPARENT_SPANGLES))&     ((self.data.cos_luz*self.data.cos_obs)<=0)
    )
    self.data.loc[cond,"transmit"]=True
    #"""
    
Spangler.update_simple_state=update_simple_state


# ### Plot observer

def plot_obs(self,show_hidden=True,center_at=None,not_plot=[],**args):
    """
    Plot spangle.

    Optional parameters:
    
        show_hidden: boolean, default = True:
            If True show hidden spangles (used to create convex hull of objects).
            
        center_at: string, default = None:
            Hash of the object around which the plotting will be centered at (see sphash column
            of the Spangler DataFrame).
            
        args: dictionary, default = dict(c='c',s=0.1):
            Scatter plotting options, dictionary.
    """
    sargs=dict(c='c',sizes=3.5)
    sargs.update(args)
    bgcolor='k'

    #Center of plot
    cond=(self.data.sphash==center_at)
    x_cen,y_cen,z_cen=self.data[cond][["x_obs","y_obs","z_obs"]].mean() if sum(cond)>0 else np.array([0,0,0])

    #Maxval original
    maxval_full=1.2*np.abs(self.data[["x_obs","y_obs","z_obs"]].to_numpy()-[x_cen,y_cen,z_cen]).max()

    #Select plotting bodies
    yes_plot=(~self.data.sphash.isin(not_plot))
    nyes_plot=sum(yes_plot)
    if nyes_plot==0:
        raise AssertionError(f"No body remain after removing {not_plot}")
    data=self.data[yes_plot]
    
    #Select scale for plot
    cond=cond if sum(cond)>0 else [True]*nyes_plot        
    maxval=1.2*np.abs(data[cond][["x_obs","y_obs","z_obs"]].to_numpy()-[x_cen,y_cen,z_cen]).max()
    size_factor=maxval_full/maxval
        
    #Figure
    fig=plt.figure(figsize=(5,5))
    fig.patch.set_facecolor(bgcolor)
    ax=fig.add_subplot(111,facecolor=bgcolor)
    ax.axis("off")

    #Plot according to state
    colors=np.array(['#000000']*nyes_plot) #Black
    sizes=np.array([0.0]*nyes_plot)

    #No illuminated
    cond=(data.visible)&(~data.illuminated)
    colors[cond]=Misc.rgb(SHADOW_COLOR_OBS,to_hex=True) #Gray
    sizes[cond]=3.5*size_factor*data.scale[cond]

    #Illuminated
    cond=(data.visible)&(data.illuminated)
    colors[cond]=[Misc.rgb([SPANGLE_COLORS[stype][0],
                            SPANGLE_COLORS[stype][1]*min((cos_luz*cos_obs+0.3),1),
                            SPANGLE_COLORS[stype][2]],
                           to_hex=True) for stype,cos_luz,cos_obs in zip(data[cond].spangle_type,
                                                                       abs(data[cond].cos_luz),
                                                                       abs(data[cond].cos_obs))
                 ] #Object color
    sizes[cond]=3.5*size_factor*data.scale[cond]

    #Transmit
    cond=(data.visible)&(data.transmit)
    colors[cond]=[Misc.rgb([SPANGLE_COLORS[stype][0],
                            SPANGLE_COLORS[stype][1]*min((cos_luz*cos_obs+0.3),1),
                            SPANGLE_COLORS[stype][2]],
                            to_hex=True) for stype,cos_luz,cos_obs in zip(data[cond].spangle_type,
                                                                       abs(data[cond].cos_luz),
                                                                       abs(data[cond].cos_obs))
                 ] #Object color
    sizes[cond]=0.5*size_factor*data.scale[cond]

    #No illuminated
    cond=(data.unset)
    colors[cond]=Misc.rgb(SHADOW_COLOR_OBS,to_hex=True) #Gray
    sizes[cond]=3.5*size_factor*data.scale[cond]
    
    #Plots
    sargs.update(dict(c=colors,sizes=sizes))    
    ax.scatter(data.x_obs-x_cen,data.y_obs-y_cen,**sargs)
    
    #Show hidden spangles
    if show_hidden:
        cond=(self.data.hidden)
        sargs.update(dict(c='r',sizes=1.5*data.scale[cond],ec='r',fc='r'))    
        ax.scatter(data.x_obs[cond]-x_cen,data.y_obs[cond]-y_cen,**sargs)

    #Ranges
    ax.set_xlim(-maxval,maxval)
    ax.set_ylim(-maxval,maxval)
    
    factor=1
    xmin,xmax=factor*np.array(list(ax.get_xlim()))
    ymin,ymax=factor*np.array(list(ax.get_ylim()))

    #Axis
    ax.plot([xmin,xmax],[0,0],'w-',alpha=0.3)
    ax.plot([0,0],[ymin,ymax],'w-',alpha=0.3)
    ax.text(xmax,0,r"$x_{obs}$",color='w',alpha=0.5,fontsize=8)
    ax.text(0,ymax,r"$y_{obs}$",color='w',alpha=0.5,fontsize=8)

    #Title
    lamb_obs=self.rqf_obs[1]*Consts.rad
    phi_obs=self.rqf_obs[2]*Consts.rad        
    label_obs=f"Obs ($\lambda$,$\\beta$) : ({lamb_obs:.1f}$^\circ$,{phi_obs:.1f}$^\circ$)"
    ax.set_title(f"Spangler {self.geometry}, N = {self.nspangles}, {label_obs}",
                 color='w',fontsize=10,position=(0.5,+0.5),ha='center')
    Plot.pryngles_mark(ax)

    #Scale
    center_text=""
    if center_at:
        center_text=f", Center at '{center_at}'"
    ax.text(0,0,f"Axis scale: {maxval*factor:.2g}{center_text}",
              fontsize=8,color='w',
              transform=ax.transAxes)

    #Decoration
    ax.axis("equal")
    fig.tight_layout()
    self.fig2d=fig
    self.ax2d=ax

Spangler.plot_obs=plot_obs


# ### Test join



# ### Update intersection state

def update_intersection_state(self):
    """Update state of intersections
    """
    self.data.intersect=True
    
    #Check if an intersection has been computed
    if len(self.qhulls) == 0:
        raise AssertionError("You must set an intersection vantage point.")
    
    for sphash in Misc.flatten([self.sphash]):
        
        cond=(self.data.sphash==sphash)
        inhull_not_in_hole=[True]
        
        verbose(VERB_SIMPLE,f"Calculating intersections for '{sphash}'")
        for i,hull in enumerate(self.qhulls[sphash]):
            
            qhull=hull["qhull"]
            if qhull is None:
                continue
            
            htype=hull["hulltype"]
            zmin,zmed,zmax=hull["zmin"],hull["zmed"],hull["zmax"]
            
            verbose(VERB_SIMPLE,f"Hull {i+1} for '{sphash}' of type '{htype}'")

            #Evaluate conditions
            inhull=sci.points_in_hull(self.data[["x_int","y_int"]],qhull)&(~cond)
            below=np.array([False]*self.nspangles)
            
            if htype=="hidden":

                #Holes
                inhull_not_in_hole=(~inhull)
                verbose(VERB_SIMPLE,f"Points outside hidden hull for '{sphash}': {sum(inhull_not_in_hole)}")
                hull["notinhole"]=sum(inhull_not_in_hole)
                continue
                
            else:
                #Body
                verbose(VERB_SIMPLE,f"Points not in hole for '{sphash}:{htype}': {sum(inhull_not_in_hole)}")

                #Spangles to evaluate
                cond_int=(~self.data.hidden)&(self.data.sphash!=sphash)&(self.data.intersect)

                if htype=="med":
                    below=(inhull_not_in_hole)&(inhull)&(self.data[cond_int]["z_int"]<zmed)
                
                elif htype=="plane":

                    #Determine if points are below
                    cond_full=(inhull_not_in_hole)&(inhull)&(cond_int)
                    verbose(VERB_SIMPLE,"Fulfilling all conditions:",sum(cond_full))
                    
                    plane=hull["plane"]
                    below[cond_full]=[not plane.is_above(r,axis=1) for r in self.data[cond_full][["x_int","y_int","z_int"]].values]
                    
                    if hull["ns_int"][2]<0:
                        below[cond_full]=~below[cond_full]
                else:
                    raise ValueError("Type of hull '{htype}' not recognized")
            
            #Store information
            verbose(VERB_SIMPLE,f"Points in hull for '{sphash}:{htype}': {sum(inhull)}")
            verbose(VERB_SIMPLE,f"Points below '{sphash}:{htype}': {sum(below)}")
            
            #Set visibility
            self.data.loc[below,"intersect"]=False

            hull["inhull"]=sum(inhull)
            hull["below"]=sum(below)
    
def update_visibility_state(self):
    self.update_intersection_state()
    self.data.visible=self.data.visible&self.data.intersect

def update_illumination_state(self):
    self.update_intersection_state()
    self.data.illuminated=self.data.illuminated&self.data.intersect

Spangler.update_intersection_state=update_intersection_state
Spangler.update_visibility_state=update_visibility_state
Spangler.update_illumination_state=update_illumination_state

def _plot_intersect(self,prop="intersect",fig=None):
    """Plot intersect
    
    Parameters:
    
        prop: string, default = intersect:
            Property to highlight.  Available: intersect, illuminated, visible.
            
        fig: Figure, default = None:
            Figure.
            
    Return:
    
        fig: Figure.
    """
    
    #Check if an intersection has been computed
    if len(self.qhulls) == 0:
        raise AssertionError("You must set an intersection vantage point.")    

    bgcolor='k'
    if fig is None:
        fig,ax=plt.subplots()
        fig.patch.set_facecolor(bgcolor)
        ax.axis("off")
    else:
        ax=fig.gca()

    #Plot qhulls and points
    for i,sphash in enumerate(Misc.flatten([self.sphash])):
        cond_obj=(self.data.sphash==sphash)
        cond=cond_obj&(~self.data.hidden)&(self.data[prop])
        
        s=1*self.data[cond].scale
        colors=[Misc.rgb(SPANGLE_COLORS[stype],to_hex=True) for stype in self.data[cond].spangle_type]
        ax.scatter(self.data[cond].x_int,self.data[cond].y_int,s=s,c=colors,zorder=100)
        
        #Plot qhull
        for hull in self.qhulls[sphash]:
            if hull["qhull"] is None:
                continue
            f=convex_hull_plot_2d(hull["qhull"],ax)

    #Remove points corresponding to qhull
    for l in fig.axes[0].get_children():
        if type(l) is Line2D:
            plt.setp(l,ms=0,zorder=100)
        if type(l) is LineCollection:
            l.set_color('w')
            
    #Decoration
    ax.set_xlabel(r"$x_{\rm int}$")
    ax.set_ylabel(r"$y_{\rm int}$")
    
    ax.grid()
    ax.axis("equal")
    
    return fig

def _view_intersect(self,lon=0,lat=0,fig=None):
    """Similar as _plot_intersect but for animation and interactive purposes.
    
    Parameters:
        long, lat: floats [deg], default = 0,0:
            Ecliptic longitud and latitude in degrees.
            
        fig: Figure, default = None.
            Figure.
    """
    
    #Input
    lon=float(lon)
    lat=float(lat)

    n_obs=sci.cartesian([1,lon*Consts.deg,lat*Consts.deg])
    cond=self.set_observer(nvec=n_obs,center=[0,0,0])

    self.reset_state()
    self.update_simple_state()
    self.update_visibility_state()
    fig=self._plot_intersect(prop="visible",fig=fig)
    
    #Decorate
    ax=fig.gca()
    ax.text(0.5,1.01,f"Observer: lon = {lon:+.2f}$^\circ$, lat = {lat:+.2f}$^\circ$",
            transform=ax.transAxes,ha='center',fontsize=10,color='w')

Spangler._plot_intersect=_plot_intersect
Spangler._view_intersect=_view_intersect

def _interact_intersect(self):
    
    verbosity=Verbose.VERBOSITY
    Verbose.VERBOSITY=VERB_NONE
    
    def view_intersect(**args):
        self._view_intersect(**args)
        
    opciones=dict(continuous_update=False,readout_format=".3f")
    interact(view_intersect,
             lon=widgets.FloatSlider(min=0,max=360,step=1,value=0,**opciones),
             lat=widgets.FloatSlider(min=-90,max=90,step=1,value=0,**opciones),
            );
    
    Verbose.VERBOSITY=verbosity

def _animate_intersect(self,filename=None,lat=20,start=0,stop=360,num=10):
    """Animate intersection state
    
    Notes:
        Based in: https://github.com/jwkvam/celluloid
    """
    
    verbosity=Verbose.VERBOSITY
    Verbose.VERBOSITY=VERB_NONE
    bgcolor='k'
    
    fig,ax=plt.subplots()
    fig.patch.set_facecolor(bgcolor)
    ax.axis("off")
    
    camera=Camera(fig)
    for lon in np.linspace(start,stop,num):
        self._view_intersect(lon=lon,lat=lat,fig=fig)
        camera.snap()
    anim=camera.animate()

    Verbose.VERBOSITY=verbosity
        
    if filename is not None:
        if 'gif' in filename:
            anim.save(filename)
            del anim
        elif 'mp4' in filename:
            ffmpeg=animation.writers["ffmpeg"]
            metadata = dict(title='Pryngles Spangler Animation',
                            artist='Matplotlib',
                            comment='Movie')
            w=ffmpeg(fps=15,metadata=metadata)
            anim.save(filename,w)
            del anim
        else:
            raise ValueError(f"Animation format '{filename}' not recognized")
    else:
        return anim

Spangler._interact_intersect=_interact_intersect
Spangler._animate_intersect=_animate_intersect


