import ctypes
import numpy as np
import glob
import pandas as pd
import time

#qtype="PLANET"
#qtype="RING_REFLECTION"
qtype="RING_TRANSMISSION"

if qtype=="PLANET":
    fileinterp="planet-interpolation.mat"
if qtype=="RING_REFLECTION":
    fileinterp="ring_forward-interpolation.mat"
if qtype=="RING_TRANSMISSION":
    fileinterp="ring_back-interpolation.mat"

#Read file with Fourier coefficients
qreflection=1
if qtype=="PLANET":
    filename="fou_gasplanet.dat"
if "RING" in qtype:
    filename="fou_ring_0_4_0_8.dat"
if "TRANSMISSION" in qtype:
    qreflection=0


#Load library
libfile = glob.glob('./cpixx*.so')[0]
cpixx=ctypes.CDLL(libfile)

DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

def double1ArrayToPointer(arr):
    arr_ptr = arr.ctypes.data_as(PDOUBLE)
    return arr_ptr
    
def double2ArrayToPointer(arr):
    """ Converts a 2D numpy to ctypes 2D array. 
    
    Arguments:
        arr: [ndarray] 2D numpy float64 array

    Return:
        arr_ptr: [ctypes double pointer]

    """

    # Init needed data types
    ARR_DIMX = DOUBLE*arr.shape[1]
    ARR_DIMY = PDOUBLE*arr.shape[0]

    # Init pointer
    arr_ptr = ARR_DIMY()

    # Fill the 2D ctypes array with values
    for i, row in enumerate(arr):
        arr_ptr[i] = ARR_DIMX()

        for j, val in enumerate(row):
            arr_ptr[i][j] = val


    return arr_ptr


def double2pointerToArray(ptr, n, m):
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

def double3ArrayToPointer(arr):
    """ Converts a 3D numpy to ctypes 3D array. 
    
    Arguments:
        arr: [ndarray] 3D numpy float64 array

    Return:
        arr_ptr: [ctypes double pointer]

    """

    # Init needed data types
    ARR_DIMX = DOUBLE*arr.shape[2]
    ARR_DIMY = PDOUBLE*arr.shape[1]
    ARR_DIMZ = PPDOUBLE*arr.shape[0]

    # Init pointer
    arr_ptr = ARR_DIMZ()

    # Fill the 2D ctypes array with values
    for i, row in enumerate(arr):
        arr_ptr[i] = ARR_DIMY()

        for j, col in enumerate(row):
            arr_ptr[i][j] = ARR_DIMX()

            for k, val in enumerate(col):
                arr_ptr[i][j][k] = val

    return arr_ptr



def double3pointerToArray(ptr, n, m, o):
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


class FourierCoefficients(ctypes.Structure):
    _fields_=[
        ("nmat",ctypes.c_int),
        ("nmugs",ctypes.c_int),
        ("nfou",ctypes.c_int),
        ("xmu",PDOUBLE),
        ("rfou",PPPDOUBLE),
        ("rtra",PPPDOUBLE),
    ]
    def __init__(self,nmat,nmugs,nfou,xmu,rfou,rtra):
        self.nmat=nmat
        self.nmugs=nmugs
        self.nfou=nfou
        self.xmu=double1ArrayToPointer(xmu)
        self.rfou=double3ArrayToPointer(rfou)
        self.rtra=double3ArrayToPointer(rtra)

#########################################
#Read fourier python
#########################################
def read_fourier(filename):
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
    f=open(filename)

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

    print(f"Checksum '{filename}': {rfou.sum()+rtra.sum():.16e}")
    return nmat,nmugs,nfou,xmu,rfou,rtra

nmat,nmugs,nfou,xmu,rfou,rtra=read_fourier(filename)
F=FourierCoefficients(nmat,nmugs,nfou,xmu,rfou,rtra)

#########################################
#Reflection, single value
#########################################
npix=1
phi=np.zeros(npix)
beta=np.zeros(npix)
theta0=np.zeros(npix)
theta=np.zeros(npix)
apix=np.zeros(npix)
Sarr=np.zeros((npix,F.nmat+1))

if qtype=="PLANET":
    phi[0]=1.0448451569439827;
    beta[0]=3.069394277348945;
    theta0[0]=0.04990329026929557;
    theta[0]=0.02509670973070432;
    apix[0]=9.432328787795567e-05;
    print("Expected values:\n4.597857424902560283e-07 2.251229972872198058e-07 1.834400800563127439e-09 4.896421313424954569e-01");
if qtype=="RING_REFLECTION":
    phi[0]=1.490116119384765625e-08;
    beta[0]=0.000000000000000000e+00;
    theta0[0]=4.999999999999998890e-01;
    theta[0]=5.000000000000000000e-01;
    apix[0]=1.163314390931110409e-04;
    print("Expected:\n6.193646058775441171e-06 -3.046132218287162406e-07 -9.759550180589223642e-15 4.918156751904256135e-02");    
if qtype=="RING_TRANSMISSION":
    phi[0]=1.601029385538801364e+00;
    beta[0]=1.601029385538801364e+00;
    theta0[0]=1.744974835125044643e-02;
    theta[0]=5.000000000000000000e-01;
    apix[0]=1.163314390931110409e-04;
    print("Expected:\n1.688771016436214060e-07 -1.485273114913191718e-08 2.674513729889046925e-09 8.936444506313186154e-02");


print(phi)
#reflection Specifications
cpixx.reflection.restype = ctypes.c_int
cpixx.reflection.argtypes = [
    FourierCoefficients,
    ctypes.c_int,
    ctypes.c_int,
    PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,
    PPDOUBLE
]
Sarr_ptr=double2ArrayToPointer(Sarr)
cpixx.reflection(F,qreflection,npix,
                 double1ArrayToPointer(phi),
                 double1ArrayToPointer(beta),
                 double1ArrayToPointer(theta0),
                 double1ArrayToPointer(theta),
                 double1ArrayToPointer(apix),
                 Sarr_ptr);
Sarr=double2pointerToArray(Sarr_ptr,*Sarr.shape)
print("Calculated values:")
print([f"{S:.16e}" for S in Sarr[0]])

#########################################
#Reflection, multiple value
#########################################
data=np.loadtxt(fileinterp)
npix=len(data)
#npix=100
phi=data[:npix,0].copy()
beta=data[:npix,1].copy()
theta0=data[:npix,2].copy()
theta=data[:npix,3].copy()
apix=data[:npix,7].copy()
Spixx=data[:npix,8:]
Sarr=np.zeros((npix,F.nmat+1))

Sarr_ptr=double2ArrayToPointer(Sarr)
cpixx.reflection(F,qreflection,npix,
                 double1ArrayToPointer(phi),
                 double1ArrayToPointer(beta),
                 double1ArrayToPointer(theta0),
                 double1ArrayToPointer(theta),
                 double1ArrayToPointer(apix),
                 Sarr_ptr);
Sarr=double2pointerToArray(Sarr_ptr,*Sarr.shape)

print("Maximum difference:",abs(Spixx-Sarr).max())
