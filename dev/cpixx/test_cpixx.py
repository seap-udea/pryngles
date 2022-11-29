import ctypes
import numpy as np
import glob
import pandas as pd
import time

libfile = glob.glob('./cpixx*.so')[0]
cpixx=ctypes.CDLL(libfile)

#########################################
#Spline
#########################################
n=10
x=np.arange(1.0,n+1,1.0)
y=x**2
y2=np.zeros_like(y)

#Specifications
#cpixx.spline.restype = np.ctypeslib.ndpointer(dtype=float)
cpixx.spline.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=float),
]
cpixx.spline(x,y,n,y2)
print(y2)

exit(0)
#########################################
#Test
#########################################
cpixx.test_cpixx()

