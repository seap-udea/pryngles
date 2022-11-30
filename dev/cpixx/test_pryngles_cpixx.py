import ctypes
import numpy as np
import glob
import pandas as pd
import time
from pryngles import *

libfile = glob.glob(Misc.get_data('../cpixx*.so'))[0]
print(libfile)

cpixx=ctypes.CDLL(libfile)
cpixx.test_cpixx()

n=10
x=np.arange(1.0,n+1,1.0)
y=x**2
y2=np.zeros_like(y)

cpixx.spline.restype = ctypes.c_double
cpixx.spline.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=float),
]

suma=cpixx.spline(x,y,n,y2)
print(y2.sum(),suma)
