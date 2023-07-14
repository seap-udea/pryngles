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
from abc import ABC, abstractmethod
from scipy.optimize import bisect
from scipy.integrate import quad,dblquad
from scipy.interpolate import interp1d,interp2d
import numpy as np
import matplotlib.pyplot as plt
import gzip

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try:
    SCATTERERS_CATALOGUE
except:
    SCATTERERS_CATALOGUE=dict()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Scatterer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Scatterer(PrynglesCommon,ABC):
    """This is an abstract class describing the basic attributes and methods
    of a Scatterer. 
    
    A Scatterer is an object able to compute the components of the Stokes vector 
    of the light scattered from a flat surface under a given ser of circumstances 
    of illumination and observation.
    
    All Scatterer has two basic properties:
    
        params: dictionary:
           A dictionary containing all the properties of the scatterer. 
    
        hash: integer: 
           An unique identification. This hash allows us to find the object in 
           a catalogue (SCATTERERS_CATALOGUE) for use at any moment in the execution.
        
    and the following general methods:
    
        register(scatterer,params):
            Register the scatterer in the Scatterers catalogue.

            Usage: always the initialization method for the daughters classes 
            must be:
            
                def __init__(self,**params):
                    if self.register(self,params):
                        ...
                        
        reset_catalogue():
            Reset the catalogue of scatterers.
            
            Usage: Scatterer.reset_catalogue()
    
    When you create a daughter class you must provide always the same method:
    
        __init__(self,**params)->int:
            Initialization of the class.

            If you use `regforce = True` the scatterer is again registered.
            
        calculate_stokes(eta:float,      -> Absolute value of the cosine of the illumination angle
                         zeta:float,     -> Absolute value of the cosine of the scattering angle
                         delta:float,    -> Difference in azimuth between incoming and observing 
                         beta:float,     -> Angle between the scattering plane and the plane of reference
                         wavelen:float,  -> Wavelength
                         reflection:bool,-> True if the light is scattered back
                         **params        -> Other parameters
                         )->float(6)   
                             
            This method calculate the stokes "albedo vector" (f,q,u,v), namely the factor by which 
            you must multiply the incoming flux Fo [W/m^2/nm or W/m^2], to obtain the stokes vector (F,Q,U,V):
            
                F = f Fo: total flux.
                Q = q Fo, U = u Fo: linearly polarized flux.
                V = v Fo: circularly polarized flux.
                
            Since normally eta, zeta, delta and beta are arrays, the routine must return a matrix 
            where each row correspond to a value of eta, etc. and the columns are the components of
            the stokes albedo vector.
            
    Example:

        You can create a Scatterer which implements this class:
    
            class MySurface(Scatterer):
            
                def __init__(self,**params):
                    if self.register(self,params):
                        #Read parameters of the scatterer
                        self.A=params["A"]
                        
                        #Initialize scatterer
                        self._initialize_scatterer()
    
                #Mandatory methods
                def calculate_stokes(self,
                                     eta,zeta,delta,beta,
                                     wavelen,reflection=True,
                                     **params):
                    f=[self.AA*eta,0,0,0]
                    return f
    
                # Private methods to prepare scatterer
                def _initialize_scatterer(self):
                    self.AA=self.A**2

    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @classmethod
    def register(self,scatterer,params):
        """Register scatterer
        """
        scatterer.params=params
        scatterer.params["name"]=scatterer.__class__.__name__
        scatterer.hash=Misc.calc_hash(params)

        #Force scatterer register
        regforce = True if 'regforce' in params else False

        if scatterer.hash in SCATTERERS_CATALOGUE:
            verbose(VERB_SIMPLE,f"Scatterer with name {scatterer.params['name']} and hash {scatterer.hash} already exist at {id(SCATTERERS_CATALOGUE)}")
            scatterer.__dict__=SCATTERERS_CATALOGUE[scatterer.hash].__dict__
            return False
        else:
            verbose(VERB_SIMPLE,f"Creating a new scatterer with name {scatterer.params['name']} and hash {scatterer.hash}")
            scatterer.params["hash"]=scatterer.hash
            SCATTERERS_CATALOGUE[scatterer.hash]=scatterer
            return True

    @classmethod
    def check_size(self,value):
        """Reset catalogue of scatterers
        """
        size = len(value) if len(np.array(value).shape)>0 else 1
        return size

    @classmethod
    def reset_catalogue(self):
        """Reset catalogue of scatterers
        """
        verbose(VERB_SIMPLE,"Resetting Scatterers catalogue")
        SCATTERERS_CATALOGUE=dict()
        
    @classmethod
    def preview_scattering(self,scatterer):

        # Ranges
        etas = np.linspace(0,1,100) 
        zetas = np.ones_like(etas) # Light is above surface
        deltas = np.zeros_like(etas)
        betas = np.zeros_like(etas)
        wavelen = 0

        # Compute reflected and transmission albedos
        fref = scatterer.calculate_stokes(etas,zetas,deltas,betas,wavelen,1)
        ftra = scatterer.calculate_stokes(etas,zetas,deltas,betas,wavelen,0)
        
        # Angles
        thetas = np.arccos(etas)*Consts.rad
        
        #Figure
        fig,axs=plt.subplots(3,2,sharex=True,figsize=(6,6))

        ax=axs[0,0]
        ax.plot(thetas,fref[:,0],'b')
            
        ax.set_ylabel("Total")
        ax.set_title("Reflection")
        
        ax=axs[1,0]
        fpol = (fref[:,1]**2+fref[:,2]**2)**0.5
        ax.plot(thetas,fpol*100,'b')
        ax.set_ylabel(fr"Polarized [$\times$ 100]")

        ax=axs[2,0]
        P = fpol/fref[:,0] if fref[:,0].sum()>0 else fpol
        ax.plot(thetas,P,'b')
        ax.set_ylabel("Degree of Polarization")
        ax.set_xlabel("Illumination angle [deg]")

        ax=axs[0,1]
        ax.plot(thetas,ftra[:,0],'r')
        ax.set_title("Transmission")

        ax.text(0.98,0.98,f"{scatterer.hash}",color='k',fontsize=8,
                ha='right',va='top',transform=ax.transAxes)

        ax=axs[1,1]
        fpol = (ftra[:,1]**2+ftra[:,2]**2)**0.5
        ax.plot(thetas,fpol*100,'r')

        ax=axs[2,1]
        P = fpol/ftra[:,0] if ftra[:,0].sum()>0 else fpol
        ax.plot(thetas,P,'r')
        ax.set_xlabel("Illumination angle [deg]")
        
        #Decoration
        for ax in Misc.flatten(axs):
            ax.grid()
            
        fig.tight_layout()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Methods to overwrite
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    @abstractmethod
    def __init__(self,**params)->str:
        pass
    
    @abstractmethod
    def calculate_stokes(self,eta:float,zeta:float,delta:float,beta:float,
                         wavelen:float,reflection:bool,
                         **params)->float:
        pass

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class NeutralSurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class WhiteSurface(Scatterer):
    """Neutral surface.
    """
    def __init__(self,**params):
        if self.register(self,params):
            pass
    
    def calculate_stokes(self,eta,zeta,delta,beta,wavelen,
                         reflection=True,**params):
        size = self.check_size(eta)
        return np.array([[1,0,0,0]]*size)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class BlackBodySurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class BlackSurface(Scatterer):
    """Black body surface
    """
    def __init__(self,**params):
        if self.register(self,params):
            pass
    
    def calculate_stokes(self,eta,zeta,delta,beta,wavelen,
                         reflection=True,**params):
        size = self.check_size(eta)
        return np.array([[0,0,0,0]]*size)
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class GraySurface
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class GraySurface(Scatterer):
    """Lambertian Gray Surface.
    
    This is the scatterer corresponding to a surface having a gray lambertian Albedo.
    
    Parameters:
    
        phase_law: function, default=lambda eta,zeta,delta,lambda:eta :

            Law of reflection (by default is Lambertian, see Russel, 1916)

            The phase_law must obey the following prototype:

                phase_law(eta,zeta,delta,**params):
                    '''Phase law of the surface

                    Parameters:
                        eta: float:
                            cosine of the incoming angle.

                        zeta: float:
                            cosine of the outgoing angle.

                        delta: float:
                            difference between the incoming and outgoing azimuth.

                        parameters: dictionary: 
                            Other parameters of the phase law.

                    Return:
                        Phase law.
                    '''
                    ...

                Other law is the Lommel-Seeliger law:

                    phase_law = lambda eta,zeta,delta:eta*zeta/(eta+zeta) (see Russell, 1916)

    """
    def __init__(self,**params):
 
        if self.register(self,params):
            verbose(VERB_SIMPLE,f"Initializing {self.params['name']} with hash {self.hash}")
            
            #Phase law
            if "phase_law" in params:
                self.phase_law=params["phase_law"]
            else:
                #Lambertian phase law
                self.phase_law=lambda eta,zeta,delta,lamb,params:eta

            #Gray albedo
            self.AL=params["AL"]

            #Calculate the gammap parameter
            self.gammap0=self._find_gammap()

            #Accelerate the calculation of the albedo
            self._accelerate_lambertian_albedo()

    def calculate_stokes(self,eta,zeta,delta,beta,wavelen,
                         reflection=True,**params):
        size = self.check_size(eta)
        f = np.zeros((size,4))

        if reflection:
            albedo = self._get_albedo(eta)
        else:
            albedo = np.zeros_like(eta)
            
        f[:,0] = albedo
        return f
    
    #####################################
    #Complimentary routines
    #####################################
    def _calc_lambertian_albedo(self,eta,gammap0=1):
        if eta==0:return self.AL
        integrand=lambda zeta:self.phase_law(eta,zeta,0,0,0)/eta
        AL=2*np.pi*gammap0*quad(integrand,0,1)[0]
        return AL

    def _find_gammap(self):
        function=lambda gammap0:self._calc_lambertian_albedo(1,gammap0)-self.AL
        gammap0=bisect(function,0.0,1.0,rtol=1e-3)
        return gammap0 if gammap0<=1 else 1
    
    def _accelerate_lambertian_albedo(self):
        etas=np.linspace(0.0,1.0,20)
        ALs=np.array([self._calc_lambertian_albedo(eta,gammap0=self.gammap0) for eta in etas])
        self._get_albedo=interp1d(etas,ALs)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class GrayAtmosphere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class GrayAtmosphere(Scatterer):
    """Gray Atmopshere.
    
    This is the scatterer corresponding to plane-parallel analytical gray atmosphere
    """
    
    def __init__(self,**params):

        if self.register(self,params):        
            verbose(VERB_SIMPLE,f"Initializing {self.params['name']} with hash {self.hash}")
            
            #Gray albedo
            self.AS=params["AS"]

            #Load reflection functions
            self._load_reflection_functions()

            #Calculate the gammap parameter
            self.gamma0=self._find_gamma()

            #Accelerate the calculation of the albedo
            self._accelerate_lambertian_albedo()

    def calculate_stokes(self,eta,zeta,delta,beta,wavelen,
                         reflection=True,**params):
        size = self.check_size(eta)
        f = np.zeros((size,4))

        if reflection:
            albedo = self._get_albedo(eta)
        else:
            albedo = np.zeros_like(eta)
            
        f[:,0] = albedo
        return f
        
    #####################################
    #Complimentary routines
    #####################################
    def _load_reflection_functions(self):
        """Load value of reflection fucntions.

        Update:
            fint: 2d interpolating function:
                x: eta (cosine incident angle)
                y: zeta (cosine scattering angle)

        Notes:
            Tab. (2.3) in Sobolev (1975).
        """
        data_ss=np.loadtxt(Misc.get_data("diffuse_reflection_function.data"))
        eta=data_ss[1:,0]
        gamma=data_ss[0,1:]
        f=data_ss[1:,1:]
        self.fint=interp2d(gamma,eta,f)  

    def _calc_reflection_coefficient(self,eta,zeta,gamma0=1):
        """Reflection coefficient of a semi-infinite (tau = infinity) atmosphere with (gray) 
        single scattering albedo gamma0

        Requires:
            - _loadReflectionFunctions

        Notes:
            Ec. (2.43) in Sobolev (1975).
        """
        rho0=gamma0*self.fint(gamma0,eta)[0]*self.fint(gamma0,zeta)[0]/(4*(eta+zeta))
        return rho0

    def _calc_spherical_albedo(self,gamma0):
        """
        Compute spherical albedo from single scattering albedo for a semi-infinite atmosphere.

        Parameters:
            gamma0: single scattering albedo (0<=gamma0<=1), float.

        Returns:
            AS: ratio of the energy diffusely reflected by a spherical planet (0<=AS<=1), float.

        Requires:
            - _loadReflectionFunctions

        Notes:
            Ec. (1.87) in Sobolev (1975).    
        """

        AS=4*dblquad(lambda y,x,*args:self._calc_reflection_coefficient(x,y,*args)*x*y,
                     0,1,lambda x:0,lambda x:1,epsrel=1e-2,args=(gamma0,))[0]
        return AS

    def _find_gamma(self):
        """
        Starting with a target spherical albedo AS, find the value of the single scattering albedo gamma0
        of a semi-infinite atmosphere having that Albedo.

        Returns:
            gamma0: the value of gamma0 corresponding to AS (0<=gamma0<=1), float.
        """
        if np.isclose(self.AS,1,rtol=1e-2):
            return 1
        function=lambda gamma0:self._calc_spherical_albedo(gamma0)-self.AS
        gamma0=bisect(function,0.0,1.0,rtol=1e-4)
        return gamma0 if gamma0<=1 else 1

    def _calc_lambertian_albedo(self,eta):
        """
        Notes: 
            Yanovistkii (1973)
        """
        integrand=lambda zeta:self._calc_reflection_coefficient(eta,zeta,gamma0=self.gamma0)*zeta
        AL=2*quad(integrand,0,1,epsrel=1e-3)[0]
        return AL

    def _accelerate_lambertian_albedo(self):
        etas=np.linspace(0,1,20)
        ALs=np.array([self._calc_lambertian_albedo(eta) for eta in etas])
        self._get_albedo=interp1d(etas,ALs)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class StokesScatterer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class FourierScatterer(Scatterer):
    """Scatterer defined in terms of a Fourier expansion of the transmission 
    and reflection matrices.
    """
    def __init__(self,**params):
        if self.register(self,params):
            verbose(VERB_SIMPLE,f"Initializing {self.params['name']} with hash {self.hash}")
            
            #Gray albedo
            self.foufile=params["foufile"]

            if not os.path.isfile(self.foufile):
                del SCATTERERS_CATALOGUE[self.hash]
                raise AssertionError(f"Fourier coefficients file '{self.foufile}' does not exist")
            
            #Load reflection functions
            self._read_foufile(self.foufile)

    def calculate_stokes(self,eta,zeta,delta,beta,wavelen,
                         reflection=True,**params):
        npix=self.check_size(eta)
        apix=np.array([1.0]*npix)
        
        #Return arrays
        Sarr=np.zeros((npix,self.F.nmat+1))
        Sarr_ptr=ExtensionUtil.mat2ptr(Sarr)
        
        #Call routine
        cpixx_ext.reflection(self.F,int(reflection),npix,
                             ExtensionUtil.vec2ptr(delta),
                             ExtensionUtil.vec2ptr(beta),
                             ExtensionUtil.vec2ptr(eta),
                             ExtensionUtil.vec2ptr(zeta),
                             ExtensionUtil.vec2ptr(apix),
                             Sarr_ptr);
        
        #Extract stokes vector
        stokes=ExtensionUtil.ptr2mat(Sarr_ptr,*Sarr.shape)
        return stokes
    
    def _read_foufile(self,filename):
        """
        Read a file containing fourier coefficients produced by PyMieDAP

        Parameters:

           filename: string:

        Returns:

            nmugs: int:
               Number of gaussian integration coefficients.

            nmat: int:
               Number of matrix.

            nfou: int:
               Number of coefficients.

            rfout: array (nmugs*nmat,nmugs,nfou):
               Matrix for the fourier coefficients for reflection.

            rtra: array (nmugs*nmat,nmugs,nfou): 
               Matrix for the fourier coefficients for transmission
        """
        is_gz = False
        with open(filename, 'rb') as test_f:
            is_gz = (test_f.read(2) == b'\x1f\x8b')

        if is_gz:
            f=gzip.open(filename,'rt')
        else:
            f=open(filename,'r')

        #Read header
        nmat=0
        imu=0
        for i,line in enumerate(f):
            if '#' in line:
                continue
            data=line.split()
            if len(data)<3:
                if len(data)==1:
                    if not nmat:
                        nmat=int(data[0])
                    else:
                        nmugs=int(data[0])
                        xmu=np.zeros(nmugs)
                else:
                    xmu[imu]=float(data[0])
                    imu+=1
            else:
                break

        #Get core data
        data=np.loadtxt(filename,skiprows=i)
        nfou=int(data[:,0].max())+1

        rfou=np.zeros((nmat*nmugs,nmugs,nfou))
        rtra=np.zeros((nmat*nmugs,nmugs,nfou))

        #Read fourier coefficients
        for row in data:
            m,i,j=int(row[0]),int(row[1])-1,int(row[2])-1
            ibase=i*nmat
            rfou[ibase:ibase+3,j,m]=row[3:3+nmat]
            if len(row[3:])>nmat:
                rtra[ibase:ibase+3,j,m]=row[3+nmat:3+2*nmat]

        verbose(VERB_SIMPLE,f"Checksum '{filename}': {rfou.sum()+rtra.sum():.16e}")
        f.close()
        
        self.nmat,self.nmugs,self.nfou=nmat,nmugs,nfou
        self.xmu,self.rfou,self.rtra=xmu,rfou,rtra
        self.F=FourierCoefficients(nmat,nmugs,nfou,xmu,rfou,rtra)
