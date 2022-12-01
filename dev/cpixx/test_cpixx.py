import ctypes
import numpy as np
import glob
import pandas as pd
import time

qtype="PLANET"
#qtype="RING_REFLECTION"
#qtype="RING_TRANSMISSION"

#Global variables
MAX_MAT=3
MAX_MUS=22
MAX_FOU=350

#Load library
libfile = glob.glob('./cpixx*.so')[0]
cpixx=ctypes.CDLL(libfile)

#Read file with Fourier coefficients
if qtype=="PLANET":
    filename="fou_gasplanet.dat"
if "RING" in qtype:
    filename="fou_ring_0_4_0_8.dat"

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

    rfou=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))
    rtra=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))

    #Read fourier coefficients
    for row in data:
        m,i,j=int(row[0]),int(row[1])-1,int(row[2])-1
        ibase=i*nmat
        rfou[ibase:ibase+3,j,m]=row[3:3+nmat]
        if len(row[3:])>nmat:
            rtra[ibase:ibase+3,j,m]=row[3+nmat:3+2*nmat]
            
    return nmat,nmugs,nfou,xmu,rfou,rtra

nmat,nmugs,nfou,xmu,rfou,rtra=read_fourier(filename)
print(f"Checksum '{filename}': {rfou.sum()+rtra.sum():.16e}")

#########################################
#Reflection
#########################################
MAX_PIX=1
phi=np.zeros(MAX_PIX)
beta=np.zeros(MAX_PIX)
theta0=np.zeros(MAX_PIX)
theta=np.zeros(MAX_PIX)
apix=np.zeros(MAX_PIX)
Sarr=np.zeros((MAX_PIX,MAX_MAT+1))
shape=np.array([nmat,nmugs,nfou],dtype=int)

npix=1
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
    
#reflection Specifications
cpixx.reflection.argtypes = [
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=int),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
]
if "TRANSMISSION" in qtype:
    cpixx.reflection(npix,phi,beta,theta0,theta,shape,xmu,rtra,apix,Sarr);
else:
    cpixx.reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);

print("Calculated values:")
print([f"{S:.16e}" for S in Sarr[0]])
