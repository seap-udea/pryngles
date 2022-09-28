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

# # Pryngles module: Science

from pryngles import *

# ## External modules

import numpy as np
import math as mh
import spiceypy as spy
from scipy.integrate import quad
from scipy.spatial import ConvexHull
from celluloid import Camera # getting the camera
import rebound as rb

# ## The Science class
# 
# The Science class is a class with routines intended to perform a wide diversity of mathematical, physical and astronomical calculations.

class Science(PrynglesCommon):pass

# ### Template method

def template(foo=1):
    """
    Method

    Parameters:

        foo: type [units], default = 1:
            Description.

    Return:

        fee: type [units]:
            Desctiption.

    """
    return foo

Science.template=template


# ### Cartesian to spherical

def spherical(xyz):
    """
    Transform cartesian coordinates into spherical coordinates

    Parameters:

        xyz: array (3):
            Cartesian coordinates

    Return:

        rqf: array (3):
            Spherical coordinates (r, theta, phi) where theta is azimutal angle and phi is 
            elevation (complement of polar angle).                

            Notice that this convention is different than that of regular vectorial calculus
            where spherical coordinates are (r,theta,phi), but theta is the polar angle and phi 
            the ezimutal one.

    """
    r,theta,phi=spy.reclat(np.array(xyz))
    theta=2*mh.pi+theta if theta<0 else theta

    return np.array([r,theta,phi])

def cospherical(xyz):
    """Transform cartesian coordinates into cosine/sine of spherical angles
    
    Parameters:

        xyz: array (3):
            Cartesian coordinates
            
    Return:

        cqsqcf: array (3):
            Cosine/sine of spherical angles (cos theta, sin theta, cos phi) where theta is 
            azimutal angle and phi is elevation (complement of polar angle).         
    """
    rho=(xyz[0]**2+xyz[1]**2)**0.5
    sf=xyz[2]/(rho**2+xyz[2]**2)**0.5
    cq=xyz[0]/rho if not mh.isclose(rho,0) else 1
    sq=xyz[1]/rho if not mh.isclose(rho,0) else 0
    return np.array([cq,sq,sf])

def pcylindrical(xyz):
    """Transform cartesian coordinates into pseudo cylindrical coordinates
    
    Parameters:

        xyz: array (3):
            Cartesian coordinates
            
    Return:

        rhoazcf: array (3):
            Cylindrical coordinates expresed as rho, phi (azimutal angle) and cos(theta) (cosine
            of polar angle).    
    """
    rho=(xyz[0]**2+xyz[1]**2)**0.5
    r=(xyz[2]**2+rho**2)**0.5
    phi=mh.atan2(xyz[1],xyz[0])
    phi=phi if phi>0 else 2*np.pi+phi
    cost=xyz[2]/r if not mh.isclose(r,0) else mh.copysign(1,xyz[2])
    return np.array([rho,phi,cost])

def cartesian(rqf):
    """
    Transform cartesian coordinates into spherical coordinates

    Parameters:

        xyz: array (3):
            Cartesian coordinates

    Return:

        rqf: array (3):
            Spherical coordinates (r, theta, phi) where theta is azimutal angle and phi is 
            elevation (complement of polar angle).                

            Notice that this convention is different than that of regular vectorial calculus
            where spherical coordinates are (r,theta,phi), but theta is the polar angle and phi 
            the ezimutal one.

    """
    return spy.latrec(rqf[0],rqf[1],rqf[2])

def direction(n):
    """Calculate the direction on which a vector is pointing
    
    Parameters:
        n: array:
            If len(n)==2, components are longitude and latitude of the direction (in degrees).
            If len(n)==3, components are cartisian coordinates of the vector n.
            
    Return:

        If len(n)==2:
        
            nvec: array(3):
                Cartesian components of the vector.

        If len(n)==3:
    
            lamb: float [degrees]:
                Longitude (angle with respect to x-axis).

            beta: float [degrees]:
                Latitude (elevation angle with respect to xy-plane).
    """
    if len(n)==3:
        rqf=spherical(n)
        return rqf[1]*Consts.rad,rqf[2]*Consts.rad
    elif len(n)==2:
        if abs(n[1])>90:
            raise ValueError("Elevation angle should be in the interval [-90,90]")
        nvec=cartesian([1,n[0]*Consts.deg,n[1]*Consts.deg])
        return nvec

Science.spherical=spherical
Science.cospherical=cospherical
Science.pcylindrical=pcylindrical
Science.cartesian=cartesian
Science.direction=direction


def rotation_matrix(ez,alpha):
    """
    Set a rotation matrix from the direction of the ez vector and a rotation angle alpha
    
    Parameter:
        ez: array (3)
            vector in the direction of the z-axis. 
            
        alpha: float (3) [rad]
            Rotation angle of the x-axis around z-axis (clockwise)
            
    Return:
        Msys2uni: float (3x3)
            Rotation matrix from the system defined by ez and the universal system.
            
        Muni2sys: float (3x3)
            Rotation matrix from the universal system to the system defined by ez
    """
    ez,one=spy.unorm(ez)
    ex=spy.ucrss([0,0,1],ez) #Spice is 5 faster for vcrss
    if spy.vnorm(ex)==0:
        ex=np.array([1,0,0]) if np.sum(ez)>0 else np.array([-1,0,0])
    ey=spy.ucrss(ez,ex)
    Msys2uni=np.array(list(np.vstack((ex,ey,ez)).transpose())).reshape((3,3))
    Muni2sys=spy.invert(Msys2uni)
    verbose(VERB_VERIFY,"Rotation axis:",ex,ey,ez)
    return Msys2uni,Muni2sys

Science.rotation_matrix=rotation_matrix


LIMB_NORMALIZATIONS=dict()

def limb_darkening(rho,cs=[0.6562],N=None):
    """
    Parameters:
        rho: float:
            Distance to center of the star in units of stellar radius.
            
        cs: list, default = [0.6562]:
            List of limb darkening coefficients.
            
        N: float, default = 1:
            Normalization constant.
            
    Return:
        I: float:
            Normalized intensity of the star at rho.
    
    Notes: 
        Models in: https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html
        Coefficients available at: https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt

    Test code:
    
        fig=plt.figure()
        ax=fig.gca()
        rhos=np.linspace(0,1,100)
        Rs=1
        coefs=[0.6550]
        N=Util.limbDarkeningNormalization(coefs)
        ax.plot(rhos,Util.limbDarkening(rhos,Rs,coefs,N))
        coefs=[0.6022,0.0654]
        N=Util.limbDarkeningNormalization(coefs)
        ax.plot(rhos,Util.limbDarkening(rhos,Rs,coefs,N))
        coefs=[0.9724,-0.4962,0.2029]
        N=Util.limbDarkeningNormalization(coefs)
        ax.plot(rhos,Util.limbDarkening(rhos,Rs,coefs,N))        
    """
    mu=(1-rho**2)**0.5
    order=len(cs)
    
    #Calculate normalization constant
    if N is None:
        chash=hash(tuple(cs))
        if chash in LIMB_NORMALIZATIONS:
            N=LIMB_NORMALIZATIONS[chash]
        else:
            integrand=lambda rho:Science.limb_darkening(rho,cs,N=1)*2*np.pi*rho
            N=quad(integrand,0.0,1.0,epsrel=1e-5)[0]
            verbose(VERB_VERIFY,f"Normalization of limb darkening function for cs = {cs}, N = {N}")
            LIMB_NORMALIZATIONS[chash]=N
            
    if order==0:
        I=np.ones_like(rho)
    elif order==1:
        I=1-cs[0]*(1-mu)
    elif order==2:
        I=1-cs[0]*(1-mu)-cs[1]*(1-mu)**2
    elif order==3:
        I=1-cs[0]*(1-mu)-cs[1]*(1-mu**1.5)-cs[2]*(1-mu**2)
    elif order==4:
        I=1-cs[0]*(1-mu**0.5)-cs[1]*(1-mu)-cs[2]*(1-mu**1.5)-cs[3]*(1-mu**2)
    else:
        raise ValueError(f"Limb darkening not implemented for order {order}")
    return I/N

Science.limb_darkening=limb_darkening


def get_convexhull(data):
    if len(data)>0:
        try:
            qhull=ConvexHull(data)
        except:
            qhull=None
    else:
        qhull=None
    return qhull
    
Science.get_convexhull=get_convexhull

def points_in_hull(p, hull, tol=1e-12):
    """Determine if a set of points are inside a convex hull.
    
    Parameters:
        
        p: numpy array (Nx2):
            Set of coordinates for points to evaluate.
    
        hull: ConvexHull:
            Convex hull to evaluate.
            
    Return:
    
        inside: boolean array (N):
            Boolean array telling if points are inside the convex hull.
            
    Examples:
        
        import numpy as np
        
        rng = np.random.default_rng()
        points = rng.random((30, 2))
        hull = ConvexHull(points)
        
        ps = rng.random((30, 2))-0.5
        cond=points_in_hull(ps,hull)

        import matplotlib.pyplot as plt
        
        for simplex in hull.simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

        for p in ps[cond]:
            plt.plot(p[0],p[1],'r*')

        for p in ps[~cond]:
            plt.plot(p[0],p[1],'co')
            
            
    Notes:
        Taken from https://stackoverflow.com/a/72483841
    """
    return np.all(hull.equations[:,:-1] @ p.T + np.repeat(hull.equations[:,-1][None,:], len(p), axis=0).T <= tol, 0)

Science.points_in_hull=points_in_hull


# ## Plane

class Plane(PrynglesCommon):
    """Get plane coefficients and coordinates.
    
    Initialization parameters:
        
        p1,p2,p3: array(3):
            Points on the plane.
            
    Attributes:
        a,b,c,d: float:
            Coefficients of the equation of the plane.
            a x + b y + c z + d = 0
            
    Notes:
        This class has been optimized removing all vectorial
        operations. This reduce considerably the execution time.
    """
    def __init__(self,p1,p2,p3):
        x1,y1,z1=p1
        x2,y2,z2=p2
        x3,y3,z3=p3
        
        a1 = x2 - x1
        b1 = y2 - y1
        c1 = z2 - z1
        a2 = x3 - x1
        b2 = y3 - y1
        c2 = z3 - z1

        self.a = b1 * c2 - b2 * c1
        self.b = a2 * c1 - a1 * c2
        self.c = a1 * b2 - b1 * a2
        self.d = (- self.a * x1 - self.b * y1 - self.c * z1)
                
        #Save components of the defining points
        self.p1x = p1[0];self.p1y = p1[1];self.p1z = p1[2]
        self.p2x = p2[0];self.p2y = p2[1];self.p2z = p2[2]
        self.p3x = p3[0];self.p3y = p3[1];self.p3z = p3[2]

        #Normal vector
        self.normal=(self.a**2+self.b**2+self.c**2)**0.5
        self.nx=self.a/self.normal
        self.ny=self.b/self.normal
        self.nz=self.c/self.normal
    
    def get_projection(self,p):
        """Find the projected point on the surface of the plane.
        
        Parameters:
            p: list (3):
                Coordinates of the point.
        
        Return:
            v: list (3):
                Coordinates of projection point.
                
            d: float:
                Distance.
        """
        
        #Distance
        d=abs(self.a*p[0]+self.b*p[1]+self.c*p[2]+self.d)/self.normal
        
        #Vectorial equivalent np.dot(p-self.p1,self.n)
        pdn=(p[0]-self.p1x)*self.nx+(p[1]-self.p1y)*self.ny+(p[2]-self.p1z)*self.nz

        #Vectorial equivalent v=p-np.dot(p-self.p1,self.n)*self.n
        v=[0]*3
        v[0]=p[0]-pdn*self.nx
        v[1]=p[1]-pdn*self.ny
        v[2]=p[2]-pdn*self.nz
        
        return v,d
    
    def get_z(self,x,y):
        """Get z value of a plane corresponding to given x, y coordinates.
        """
        z = (-self.a*x-self.b*y-self.d)/self.c if not mh.isclose(self.c,0) else np.nan
        return z
    
    def is_above(self,p,vdir):
        """Check if a point is above or below a plane with respect to a given direction
        
        Parameters:
            p: list (3):
                Coordinates of the point.
                
            vidr: list (3):
                Direction with respect to
        """
        v,d=self.get_projection(p)        
        #Sign of (v-p).vdir
        cdir=(v[0]-p[0])*vdir[0]+(v[1]-p[1])*vdir[1]+(v[2]-p[2])*vdir[2]
        return cdir<=0
    
    def is_below(self,p,vdir):
        return not self.is_above(p,vdir)
    
    def plot_plane(self,ax=None,p=None,**args):
        
        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(111,projection='3d')
        
        maxval=max(abs(self.p1x),abs(self.p1y),abs(self.p1z),
                   abs(self.p2x),abs(self.p2y),abs(self.p2z),
                   abs(self.p3x),abs(self.p3y),abs(self.p3z))
        
        if p is not None:
            maxval=max(maxval,abs(p[0]),abs(p[1]),abs(p[2]))

        X,Y = np.meshgrid(np.linspace(-maxval,+maxval),np.linspace(-maxval,+maxval))
        Z=self.get_z(X,Y)
        
        ax.plot_surface(X,Y,Z,**args)
        ax.plot([self.p1x,self.p2x,self.p3x],
                [self.p1y,self.p2y,self.p3y],
                [self.p1z,self.p2z,self.p3z],'co')
        
        f=2
        ax.set_xlim(-f*maxval,+f*maxval)
        ax.set_ylim(-f*maxval,+f*maxval)
        ax.set_zlim(-f*maxval,+f*maxval)
        
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        if ax is None:
            ax.set_title(f"Plane: $a = {self.a:.2f}$, $b = {self.b:.2f}$, $c = {self.c:.2f}$, $d = {self.d:.2f}$")
        
        self.ax=ax
        self.maxval=maxval
        
        if p is not None:
            self.plot_point(p)
        
    def plot_point(self,p):
        maxval=max(self.maxval,abs(p[0]),abs(p[1]),abs(p[2]))
        ax=self.ax
        
        ax.plot(p[0],p[1],p[2],'ro')
        v,d=self.get_projection(p)
        ax.plot(v[0],v[1],v[2],'b*')

        ax.plot([p[0],v[0]],
                [p[1],v[1]],
                [p[2],v[2]],
                'r--')

        f=2
        ax.set_xlim(-f*maxval,+f*maxval)
        ax.set_ylim(-f*maxval,+f*maxval)
        ax.set_zlim(-f*maxval,+f*maxval)

Science.Plane=Plane


# ## Flybys

def calc_flyby(normal=[0,0,1],start=0,stop=360,num=10,lat=0):
    
    """Calculate a flyby coordinates
    
    Parameters:
        normal: array (3), default = [0,0,1]:
            Normal to flyby plane.
            
        start: float, default = 0:
            Start longitude.
            
        stop: float, default = 0:
            Stop longitude.
            
        num: int, default = 10:
            Number of points in flyby.
            
        lat: float, default = 0:
            Constant latitude of flyby.
    """

    #Range of longitudes and latitudes
    lonp=np.linspace(start,stop,num)
    latp=lat*np.ones_like(lonp)
    
    #Rotation matrices
    M,I=Science.rotation_matrix(normal,0)

    #Compute directions
    nvecs=np.zeros((num,3))
    for i in range(num):
        rp=Science.direction([lonp[i],latp[i]])
        nvecs[i]=spy.mxv(I,rp)

    return nvecs

Science.calc_flyby=calc_flyby


# ## Animate rebound

def animate_rebound(sim,tini=0,tend=2*np.pi,nsnap=10,interval=100):

    verbosity=Verbose.VERBOSITY
    Verbose.VERBOSITY=VERB_NONE
    
    fig,ax=plt.subplots()

    camera=Camera(fig)

    for t in np.linspace(tini,tend,nsnap):
        sim.integrate(t)
        sim.move_to_com()
        
        for p in sim.particles:
            xyz=p.xyz
            ax.plot(xyz[0],xyz[1],marker='o',color='r')
        
        camera.snap()
    
    ax.axis("equal")
    ax.grid()
    
    anim=camera.animate(interval=interval)
    
    Verbose.VERBOSITY=verbosity
    
    return anim

Science.animate_rebound=animate_rebound


