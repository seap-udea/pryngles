import ctypes
import numpy as np
import glob
import pandas as pd
import time

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

#Real test


#Test routines
n=3
m=3
p=3
xmu=np.random.rand(m)
rfou=np.random.rand(n,m,p)
rtra=np.random.rand(n,m,p)
print(rfou)
F=FourierCoefficients(n,m,p,xmu,rfou,rtra)
cpixx.sum_structure.restype = ctypes.c_double
cpixx.sum_structure.argtypes = [FourierCoefficients,ctypes.c_int,ctypes.c_int,ctypes.c_int]
suma=cpixx.sum_structure(F,n,m,p)
print(suma)
C=double3pointerToArray(F.rfou,*rfou.shape)
print(C)
exit(0)

n=3
m=3
p=3
C=np.random.rand(n,m,p)
Cp=double3ArrayToPointer(C)
print(C)
cpixx.sum_cube.restype = ctypes.c_double
cpixx.sum_cube.argtypes = [PPPDOUBLE,ctypes.c_int,ctypes.c_int,ctypes.c_int]
suma=cpixx.sum_cube(Cp,n,m,p)
print(f"Suma = {suma}");
C=double3pointerToArray(Cp,*C.shape)
print(C)
exit(0)

n=3
m=3
#M=np.ones((n,m))
M=np.random.rand(n,m)
print(M)
Mp=double2ArrayToPointer(M)
cpixx.sum_matrix.restype = ctypes.c_double
cpixx.sum_matrix.argtypes = [PPDOUBLE,ctypes.c_int,ctypes.c_int]
suma=cpixx.sum_matrix(Mp,n,m)
M=double2pointerToArray(Mp,*M.shape)
print(f"Suma = {suma}");
print(M)
exit(0)

n=30
v=np.arange(1.0,n,1.0)
print(v)
vp=double1ArrayToPointer(v)
cpixx.sum_vector.restype = ctypes.c_double
cpixx.sum_vector.argtypes = [PDOUBLE,ctypes.c_int]
suma=cpixx.sum_vector(vp,n)
print(f"Suma = {suma}");
print(v)
exit(0)

