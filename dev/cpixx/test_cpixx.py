import ctypes
import numpy as np
import glob
import pandas as pd
import time

libfile = glob.glob('./cpixx_dev*.so')[0]
cpixx=ctypes.CDLL(libfile)

filename=b"fou_gasplanet.dat";nmat=3;nmugs=21;nfou=3
#filename=b"fou_ring_0_4_0_8.dat";nmat=3;nmugs=20;nfou=344
MAX_MAT=3
MAX_MUS=22
MAX_FOU=350

#########################################
#Dynamic allocation of matrix
#########################################
M=np.ones((3,3))
suma=cpixx.sum_matrix(M,3,3);
exit(0);

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
    nfou=int(data[:,0].max())

    #rfou=np.zeros((nmat*nmugs,nmugs,nfou+1))
    #rtra=np.zeros((nmat*nmugs,nmugs,nfou+1))
    rfou=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))
    rtra=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))

    #Read fourier coefficients
    for row in data:
        m,i,j=int(row[0]),int(row[1])-1,int(row[2])-1
        ibase=i*nmat
        rfou[ibase:ibase+3,j,m]=row[3:3+nmat]
        if len(row[3:])>nmat:
            rtra[ibase:ibase+3,j,m]=row[3+nmat:3+2*nmat]
            
    return nmat,nmugs,xmu,rfou,rtra

nmat,nmugs,xmu,rfou,rtra=read_fourier(filename.decode())
print(nmat,nmugs,xmu)
print(f"Checksum rfou: {rfou.sum()}")
print(f"Checksum rtra: {rtra.sum()}")

#########################################
#Read fourier
#########################################
"""
shape=np.zeros(3,dtype=int)
#xmu=np.zeros(nmugs)
#rfou=np.zeros((nmat*nmugs,nmugs,nfou))
#rtra=np.zeros((nmat*nmugs,nmugs,nfou))
xmu=np.zeros(MAX_MUS)
rfou=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))
rtra=np.zeros((MAX_MAT*MAX_MUS,MAX_MUS,MAX_FOU))

#read_fourier Specifications
cpixx.read_fourier.argtypes = [
    ctypes.c_char_p,
    np.ctypeslib.ndpointer(dtype=int),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
]
cpixx.read_fourier(filename,shape,xmu,rfou,rtra)
print(f"Python checksum rfou: {rfou.sum():e}")
print(f"Python checksum rtra: {rtra.sum():e}")
print(shape)
for i in range(nmugs*nmat):
    for j in range(nmugs):
        for m in range(nfou):
            print(f"{i} {j} {m} {rfou[i,j,m]:.16e}")
"""

#########################################
#Reflection
#########################################
MAX_PIX=1000
phi=np.zeros(MAX_PIX)
beta=np.zeros(MAX_PIX)
theta0=np.zeros(MAX_PIX)
theta=np.zeros(MAX_PIX)
apix=np.zeros(MAX_PIX)
Sarr=np.zeros((MAX_PIX,MAX_MAT+1))
#"""
shape=np.array([nmat,nmugs,nfou],dtype=int)
#"""
print(shape)

npix=1
phi[0]=1.0448451569439827;
beta[0]=3.069394277348945;
theta0[0]=0.04990329026929557;
theta[0]=0.02509670973070432;
apix[0]=9.432328787795567e-05;

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
cpixx.reflection(npix,phi,beta,theta0,theta,shape,xmu,rfou,apix,Sarr);

print(Sarr[0])
exit(0)

#########################################
#Spline
#########################################
n=10
x=np.arange(1.0,n+1,1.0)
y=x**2
y2=np.zeros_like(y)

#Spline Specifications
#cpixx.spline.restype = np.ctypeslib.ndpointer(dtype=float)
cpixx.spline.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=float),
]
cpixx.spline(x,y,n,y2)
print(y2)

#########################################
#Splint
#########################################
#Splint Specifications
cpixx.splint.restype = ctypes.c_double
cpixx.splint.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
    ctypes.c_double,
]

for i in range(1,n):
    xv=i+0.3;
    yv=cpixx.splint(x,y,y2,n,xv);
    print(xv,yv);

#########################################
#Test
#########################################
cpixx.test_cpixx()

