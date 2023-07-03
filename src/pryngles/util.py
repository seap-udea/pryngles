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
import time
import inspect
import os
import re
import gdown
import ctypes
import glob
import math
import pandas as pd
import numpy as np
import math as mh
import spiceypy as spy
import rebound as rb
from tqdm import tqdm
from celluloid import Camera # Deprecated
from colorsys import hls_to_rgb
from collections import OrderedDict as odict
import matplotlib.pyplot as plt
from collections.abc import Iterable
from scipy.integrate import quad
from scipy.spatial import ConvexHull

#Plotting in 3d
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch
from mpl_toolkits import mplot3d
from scipy.spatial.transform import Rotation
from matplotlib import animation

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constants 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hash maxsize
from sys import maxsize as HASH_MAXSIZE

# C types for extensions
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

# Set of limb darkening normalizations
SCIENCE_LIMB_NORMALIZATIONS=dict()

# Class constants
class Consts(object):
    """Constants class
    """
    def get_physical():
        import pryngles as pr
        all_constants=[]
        for key in Consts.__dict__.keys():
            patterns = "^[a-z]+$"
            if re.search(patterns,key):
                all_constants+=[key]
        return sorted(all_constants)

    def get_all():
        import pryngles as pr
        all_constants=[]
        for key in pr.__dict__.keys():
            patterns = "^[A-Z_]+$"
            if re.search(patterns,key):
                all_constants+=[key]
        return sorted(all_constants)

#Mathematical constants
Consts.rad=180/np.pi
Consts.deg=1/Consts.rad
Consts.ppm=1e6 #parts per million factor
Consts.ppb=1e9 #parts per billion factor

#Physical constants
GSI=rb.units.convert_G(["m","s","kg"]) # G constant in SI units
for const in "times","lengths","masses":
    values=eval(f"rb.units.{const}_SI.copy()")
    for key in values:
        exec(f"Consts.{key}=values[key]")

#Size of reference objects
Consts.rearth=6378.137e3 #m, volumetric mean radius, source: 
Consts.rsun=695700e3 #m, nominal solar radius, source: 
Consts.rjupiter=71492e3 #m, equatorial radius, source: 
Consts.rsaturn=60268e3 #m, equatorial radius, source: 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Misc
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Misc(object):
    """
    Miscelaneous routines.
    
    This is a set of util routines intended for a diversity of purposes.
    
    Routines included:
    
        get_data(file)
    """

    def convert_unit(value,
                     units=dict(us=1e-6,ms=1e-3,s=1e0,m=60,h=3600)):
        """
        Convert a value to the neares multiple and submultiple
        """
        minimum=1e100
        quantity=-1
        for u,f in units.items():
            convert=value/f
            if (convert>1) and (convert<minimum):
                quantity=convert
                unit=u
                minimum=convert
        if quantity<0:
            unit=list(units.keys())[0]
            quantity=value/units[unit]
        return quantity,unit
        
    TIME_IMPORT=-1
    TIME_FIRST=-1
    TIME_LAST=-1
    def elapsed_time(show=True, total=False, imp=False, msg=None,
                     verbosity=VERB_NONE):
        """Calculate the time elapsed from last call
        Parameters:

            show: boolean, default = True:
                If true show the elapsed time.
                If false just mark the time.

            total: boolean, default = False:
                If true compute the elapsed time since the first call of the
                function

            imp: boolean, default = False:
                If true compute the elapsed time since the first call of the
                function

        """
        if Misc.TIME_IMPORT<0:
            Misc.TIME_IMPORT=time.time()     
        elif Misc.TIME_FIRST<0:
            Misc.TIME_FIRST=time.time()
            Misc.TIME_LAST=time.time()
        elif Misc.TIME_LAST<0:
            Misc.TIME_LAST=time.time()
        
        tnow=time.time()
        if total:
            tpast = Misc.TIME_FIRST
            ref = 'first call'
        elif imp:
            tpast = Misc.TIME_IMPORT
            ref = 'import'
        else:
            tpast = Misc.TIME_LAST
            Misc.TIME_LAST=tnow
            ref = 'last call'

        if show:
            value=tnow-tpast
            elapsed,unit=Misc.convert_unit(value)
            if msg is None:
                mgs = f"Elapsed time since {ref}:"
            verbose(verbosity,f"{msg}: {value:g} s = {elapsed:.2f} {unit}")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Data methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def get_data(path):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            datafile: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(os.path.dirname(__file__),'data',path);

    def retrieve_data(datafile,path="/tmp/",quiet=False,overwrite=False):
        """Retrieve a data file from public Pryngles repo, https://bit.ly/pryngles-data.

        Parameters:
        
            datafile: string.
                Name of the datafile to retrieve.

            path: string, default = '/tmp/'
                Path where the data files be retrieved.

            quiet: bool, default = False
                Is the downloading process quiet or not. This is a `gdown` option.

        Return:
            List of datafiles retrieve it.        
        """
        # Get the list of files
        filename = path+'/'+DATA_INDEX['filename']
        if not os.path.isfile(filename) or overwrite:
            url = DATA_INDEX['baseurl']+DATA_INDEX['fileid']
            gdown.download(url,filename,quiet=quiet)
        else:
            if not quiet:
                print(f"Index file {filename} already retrieved. For overwrite use overwrite = True.")
            
        files=pd.read_csv(filename,index_col=DATA_INDEX['col_filename'])
        if not quiet:
            print(f"There are {len(files)} files in data repository.")

        if isinstance(datafile,str):
            datafile=[datafile]

        if len(datafile) == 0:
            return list(files.index)
            
        dfiles=[]
        for dfile in datafile:
            # Look for file in the data index
            if dfile in files.index:
                fileid = files.loc[dfile,DATA_INDEX['col_fileid']]
                # Download
                url = DATA_INDEX['downurl']+fileid
                filename = path+'/'+dfile
                if not os.path.isfile(filename) or overwrite:
                    gdown.download(url,filename,quiet=quiet)
                else:
                    if not quiet:
                        print(f"File {filename} already retrieved. For overwrite use overwrite = True.")
                dfiles+=[filename]
            else:
                raise ValueError(f"Datafile {dfile} not available in repository. List of available files:\n{list(files.index)}")

        if not quiet:
            print(f"Files downloaded: {dfiles}")

        return dfiles

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Input/output methos
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def print_df(df):
        """Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        from IPython.display import display,HTML
        display(HTML(df.to_html()))
        
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Array methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def flatten(collection):
        """Flatten a list of objects

        Examples:
            list(Misc.flatten(["cosa"]))
            list(Misc.flatten([["cosa"]]))
            list(Misc.flatten([["cosa","perro"]]))
            list(Misc.flatten([[1,"perro"],object,float]))
        """
        for i in collection:
            if isinstance(i, Iterable) and not isinstance(i, basestring):
                for subc in Misc.flatten(i):
                    yield subc
            else:
                yield i
                
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Programming methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def get_methods(my_class):
        """Get a list of the methods for class my_class
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])
    
    def calc_hash(obj):
        if type(obj) is dict:
            hash_obj=frozenset(obj.items())
        else:
            hash_obj=obj
        hash_val=str(hash(hash_obj)%((HASH_MAXSIZE+1)*2))
        return hash_val

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class ExtensionUtil
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class ExtensionUtil(object):
    """Util routines for extensions.
    """
    def vec2ptr(arr):
        """Converts a 1D numpy to ctypes 1D array. 

        Parameters:
            arr: [ndarray] 1D numpy float64 array

        Return:
            arr_ptr: [ctypes double pointer]
        """
        arr_ptr = arr.ctypes.data_as(PDOUBLE)
        return arr_ptr

    def mat2ptr(arr):
        """ Converts a 2D numpy to ctypes 2D array. 

        Arguments:
            arr: [ndarray] 2D numpy float64 array

        Return:
            arr_ptr: [ctypes double pointer]

        """

        ARR_DIMX = DOUBLE*arr.shape[1]
        ARR_DIMY = PDOUBLE*arr.shape[0]

        arr_ptr = ARR_DIMY()

        # Fill the 2D ctypes array with values
        for i, row in enumerate(arr):
            arr_ptr[i] = ARR_DIMX()

            for j, val in enumerate(row):
                arr_ptr[i][j] = val

        return arr_ptr

    def ptr2mat(ptr, n, m):
        """ Converts ctypes 2D array into a 2D numpy array. 

        Arguments:
            arr_ptr: [ctypes double pointer]

        Return:
            arr: [ndarray] 2D numpy float64 array

        """

        arr = np.zeros(shape=(n, m))

        for i in range(n):
            for j in range(m):
                arr[i,j] = ptr[i][j]

        return arr

    def cub2ptr(arr):
        """ Converts a 3D numpy to ctypes 3D array. 

        Arguments:
            arr: [ndarray] 3D numpy float64 array

        Return:
            arr_ptr: [ctypes double pointer]

        """

        ARR_DIMX = DOUBLE*arr.shape[2]
        ARR_DIMY = PDOUBLE*arr.shape[1]
        ARR_DIMZ = PPDOUBLE*arr.shape[0]

        arr_ptr = ARR_DIMZ()

        # Fill the 2D ctypes array with values
        for i, row in enumerate(arr):
            arr_ptr[i] = ARR_DIMY()

            for j, col in enumerate(row):
                arr_ptr[i][j] = ARR_DIMX()

                for k, val in enumerate(col):
                    arr_ptr[i][j][k] = val

        return arr_ptr

    def ptr2cub(ptr, n, m, o):
        """ Converts ctypes 3D array into a 3D numpy array. 

        Arguments:
            arr_ptr: [ctypes double pointer]

        Return:
            arr: [ndarray] 3D numpy float64 array

        """

        arr = np.zeros(shape=(n, m, o))

        for i in range(n):
            for j in range(m):
                for k in range(o):
                    arr[i,j,k] = ptr[i][j][k]

        return arr

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Science
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Science(PrynglesCommon):
    """Science utility Class
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file science
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    def direction(*args):
        """Calculate the direction on which a vector is pointing
        
        Parameters:
            args: array:
                If len(args)==2, components are longitude and latitude of the direction (in degrees).
                If len(args)==3, components are cartisian coordinates of the vector n.
                
        Return:
    
            If len(args)==2:
            
                vx,vy,vz: float:
                    Cartesian components of the vector.
    
            If len(args)==2:
        
                lamb: float [degrees]:
                    Longitude (angle with respect to x-axis).
    
                beta: float [degrees]:
                    Latitude (elevation angle with respect to xy-plane).
        """
        if len(args)==3:
            rqf=Science.spherical(list(args))
            return rqf[1]*Consts.rad,rqf[2]*Consts.rad
        elif len(args)==2:
            if abs(args[1])>90:
                raise ValueError("Elevation angle should be in the interval [-90,90]")
            nvec=Science.cartesian([1,args[0]*Consts.deg,args[1]*Consts.deg])
            return nvec
        else:
            raise ValueError("You provided a wrong number of arguments '{args}'.  It should be 2 or 3'")
    
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
    
    def limb_darkening(rho,cs=[0.6562],N=None):
        """Limb darkening computation
        
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
            if chash in SCIENCE_LIMB_NORMALIZATIONS:
                N=SCIENCE_LIMB_NORMALIZATIONS[chash]
            else:
                integrand=lambda rho:Science.limb_darkening(rho,cs,N=1)*2*np.pi*rho
                N=quad(integrand,0.0,1.0,epsrel=1e-5)[0]
                verbose(VERB_VERIFY,f"Normalization of limb darkening function for cs = {cs}, N = {N}")
                SCIENCE_LIMB_NORMALIZATIONS[chash]=N
                
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
    
    def get_convexhull(data):
        if len(data)>0:
            try:
                qhull=ConvexHull(data)
            except:
                qhull=None
        else:
            qhull=None
        return qhull
        
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Plane
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Plane(PrynglesCommon):
    """A plane in 3d.
    
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Plot(object):
    """Plotting util class
    """
    
    def _pathpatch_2d_to_3d(pathpatch,pivot=[0,0,0],zDir=[0,0,1]):
        """
        Create a patch in 3d around pivot in the direction of zDir
        
        Source: https://stackoverflow.com/a/69785236
        """

        path = pathpatch.get_path() #Get the path and the associated transform
        trans = pathpatch.get_patch_transform()
        path = trans.transform_path(path) #Apply the transform

        pathpatch.__class__ =  mplot3d.art3d.PathPatch3D #Change the class
        pathpatch._path2d = path       #Copy the 2d path
        pathpatch._code3d = path.codes #Copy the codes
        pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color

        # Get the 2D vertices and add the third dimension
        verts3d = np.empty((path.vertices.shape[0],3))
        verts3d[:,0:2] = path.vertices
        verts3d[:,2] = pivot[2]

        #Get rotation matriz
        norm = np.linalg.norm(zDir)
        zDir = zDir/norm
        if np.abs(zDir[2])==1:
            yDir = np.array([0,zDir[2],0])
        else:
            yDir = (np.array([0,0,1]) - zDir[2]*zDir)/math.sqrt(1-zDir[2]**2)
        rotMat = np.empty((3,3))
        rotMat[:,0] = np.cross(zDir,yDir)
        rotMat[:,1] = yDir
        rotMat[:,2] = -zDir
        R=Rotation.from_matrix(rotMat)

        #Displace
        pathpatch._segment3d = R.apply(verts3d - pivot) + pivot

        return pathpatch

    # places a 3D circle in axes with 3d projection. 
    def circle3d(ax, center = (0,0,0), radius = 1, zDir='z', **kwargs):
        """
        Add a circle in 3d
        """
        pc = Circle(center[0:2], radius, **kwargs)
        ax.add_patch(Plot._pathpatch_2d_to_3d(pc, center, zDir))
        
    def pryngles_mark(ax):
        """Add a water mark to a 2d or 3d plot.
        
        Parameters:
        
            ax: Class axes: 
                Axe where the pryngles mark will be placed.
        """
        #Get the height of axe
        axh=ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height
        fig_factor=axh/4
        
        #Options of the water mark
        args=dict(
            rotation=270,ha='left',va='top',
            transform=ax.transAxes,color='pink',fontsize=8*fig_factor,zorder=100
        )
        
        #Text of the water mark
        mark=f"Pryngles {version}"
        
        #Choose the according to the fact it is a 2d or 3d plot
        try:
            ax.add_collection3d
            plt_text=ax.text2D
        except:
            plt_text=ax.text
            
        text=plt_text(1,1,mark,**args);
        return text

    def rgb(hls,to_hex=False):
        """Convert from hue (0-360), level (0-1) and saturation (0-1) to RGB
        
        Parameters:
        
            hls: array(3):
                Array with values of color:
                    hls[0]: hue, 0-360, see https://pythonfordesigners.com/tutorials/hsl-color-wheel/
                    hls[1]: level, 0: black, 1: white
                    hls[2]: saturation, 0: gray, 1: full-color
                    
        Return:
        
            rgb: array(3):
                Array with rgb values (R: red, G: green, B: blue)
        """
        rgb_color=hls_to_rgb(hls[0]/360.0,hls[1],hls[2])
        if to_hex:
            hex_color="#{:02x}{:02x}{:02x}".format(int(rgb_color[0]*255),
                                                   int(rgb_color[1]*255),
                                                   int(rgb_color[2]*255))
            return hex_color
        return rgb_color
    
    def rgb_sample(H=0):
        """Create a color table for a given hue
        """
        fig,ax=plt.subplots(figsize=(9,9))
        dL=0.1
        dS=0.1
        for S in np.arange(0,1+dS,dS):
            for L in np.arange(0,1+dL,dL):
                c=Circle((L,S),dL/2.5,color=Plot.rgb([H,L,S]))
                ax.add_patch(c)
                ax.text(L,S,f"S={S:.1g},L={L:.1g}",ha='center',va='center',fontsize=6,color='y')
        ax.axis("off")
        ax.axis("equal")
        plt.tight_layout()

    def draw_pryngles(iobs=1,dark=False):
        ############################################################
        # Simulation
        ############################################################
        no=10000
        P=RingedPlanet(Nr=500,Np=500,Nb=0,
                       Rint=1.2,Rext=2.0, i=45*DEG,
                       a=0.1,e=0.1,lambq=70*DEG,
                       physics=dict(AL=1,AS=1,taug=1),
                       behavior=dict(shadows=0))
        P.changeObserver([90*DEG,iobs*DEG]) #LOGO
        lamb_initial = 0.0*DEG
        lamb_final   = 360*DEG
        lambs        = np.linspace(lamb_initial,lamb_final,no)
        Rps=[]
        Rrs=[]
        ts=[]
        Ts =[]
    
        for lamb in lambs:
            P.changeStellarPosition(lamb)
            ts+=[P.t*P.CU.UT]
            P.updateOpticalFactors()
            P.updateDiffuseReflection()
            P.updateTransit()
    
            Rps+=[P.Rip.sum()]
            Rrs+=[P.Rir.sum()]
            Tp=P.Tip.sum()
    
            T=Tp+P.Tir.sum()
            Ts+=[T]
    
        Ts=np.array(Ts)
        ts=np.array(ts)
        Rps=np.array(Rps)
        Rrs=np.array(Rrs)
        ts=ts/Const.days
    
        ############################################################
        # Plot
        ############################################################
        ppm=1e6
        alpha=1
    
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(projection='3d')
        title=f"$a={P.a:g}$ au, $i={P.i*RAD:g}^\circ$ ($i_0={P.io*RAD:.1f}^\circ$), $\lambda_\mathrm{{q}}={P.lambq*RAD:g}^\circ$, Obs ($\lambda$,$\\beta$) : ({P.eobs_ecl[0]*RAD:g}$^\circ$,{P.eobs_ecl[1]*RAD:g}$^\circ$)"
        theta= np.linspace(0, 2*np.pi, no)
        x    = np.cos(theta)
        y    = np.sin(theta)
        z1    = (ppm*(Rrs+Rps-1e-3*Ts))
        z2    = (ppm*(Rps+Rps-1e-3*Ts))
    
        if dark:
            cmap=cmr.bubblegum
            back='k'
            fcolor='pink'
        else:
            cmap=cmr.bubblegum_r
            back='w'
            fcolor='blue'
    
        p=ax.scatter(x,y,(z2), marker=".", c=z1+z2 ,s=7, cmap=cmap,alpha=alpha,edgecolors=None)
        cb=plt.colorbar(p, orientation="horizontal", fraction=0.03, pad=-0.2)
        cb.set_label(r"Pryngles", color=fcolor, fontsize=40,fontname="Special Elite")
    
        cb.ax.tick_params(labelcolor=back)
        cbytick_obj = plt.getp(cb.ax, 'yticklabels' ) #Set y tick label color
        plt.setp(cbytick_obj, color=back)
        cb.ax.tick_params(which = 'minor', length = 2, color = back )
        cb.ax.tick_params(which = 'major', length = 5, color = back )
        cb.update_ticks()
    
        # THE (SPHERICAL) SUN 
        ax.scatter(0,0,0, marker=".", s=1000, color="orange")
        for i in np.linspace(0,1,20):
            ax.scatter(0,0,0, marker=".", s=1000*5*i, color="gold", alpha=1-0.9*i )
    
        # AXIS SETUP
        fig.set_facecolor(back)
        ax.set_axis_off()
        ax.set_facecolor(back) 
    
        ##CAMERA ORIENTATION 
        ax.view_init(elev=-30, azim=25)
        fig.tight_layout()
    

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
            rp=Science.direction(lonp[i],latp[i])
            nvecs[i]=spy.mxv(I,rp)
    
        return nvecs
    
    def animate_rebound(sim,filename=None,tini=0,tend=None,nsnap=None,interval=100,axis=False,traces=False,**plot_args):
        """Animate a rebound simulation.
        """
        default_plot_args=dict(
            marker='o',
            color='r'
        )
        default_plot_args.update(plot_args)
        
        verbosity=Verbose.VERBOSITY
        Verbose.VERBOSITY=VERB_NONE
        
        fig,ax=plt.subplots()
    
        if not traces:
            camera=Camera(fig)
    
        #Get the period of the longest osculant orbit
        P=-1
        for p in sim.particles[1:]:
            P=p.P if p.P>P else P
        
        #Choose properly tend and nsnap
        tend=P if tend is None else tend
        nsnap=int(tend/(P/100)) if nsnap is None else nsnap
        
        if traces:
            sim.move_to_com()
            for p in sim.particles:
                xyz=p.xyz
                ax.plot(xyz[0],xyz[1],marker="*",color='k',ms=10,zorder=1000)
    
        #Simulate
        for i,t in enumerate(tqdm(np.linspace(tini,tend,nsnap))):
            sim.integrate(t)
            sim.move_to_com()
            
            for p in sim.particles:
                xyz=p.xyz
                ax.plot(xyz[0],xyz[1],**default_plot_args)
             
            if not traces:
                ax.text(0.5,1,f"t = {sigfig.round(t,3)} (snap {i+1}/{nsnap})",
                        transform=ax.transAxes,
                        ha='center',va='bottom')
    
                camera.snap()
        
        if axis:
            ax.grid()
        else:
            ax.axis("off")
        ax.axis("equal")
    
        if not traces:
            anim=camera.animate(interval=interval)    
            Verbose.VERBOSITY=verbosity
    
            if filename is not None:
                if 'gif' in filename:
                    anim.save(filename)
                    return anim
                elif 'mp4' in filename:
                    ffmpeg=animation.writers["ffmpeg"]
                    metadata = dict(title='Pryngles Spangler Animation',
                                    artist='Matplotlib',
                                    comment='Movie')
                    w=ffmpeg(fps=15,metadata=metadata)
                    anim.save(filename,w)
                    return anim
                else:
                    raise ValueError(f"Animation format '{filename}' not recognized")
            else:
                return anim
    
            return anim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Initialization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Load library
libfile = glob.glob(Misc.get_data('../cpixx*.so'))[0]
cpixx_ext=ctypes.CDLL(libfile)

cpixx_ext.reflection.restype = ctypes.c_int
cpixx_ext.reflection.argtypes = [
    ctypes.Structure,
    ctypes.c_int,
    ctypes.c_int,
    PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,
    PPDOUBLE
]
