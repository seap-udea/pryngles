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
#Input 
n=10
x=np.arange(1.0,n+1,1.0)
y=x**2
y2=np.zeros_like(y)

data=pd.DataFrame(np.array([x,y,y2]).T,columns=["x","y","y2"])

x_arg=ctypes.c_double*n
suma=cpixx.spline(x_arg(*x),x_arg(*y),n,x_arg(*y2))
print(y2)
print(suma)

cpixx.spline.restype = ctypes.c_double
cpixx.spline.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=float),
]
suma=cpixx.spline(x,y,n,y2)
print(y2.sum(),suma)

st=time.process_time()
suma=cpixx.spline(x,y,n,y2)
data["y2"]=y2
et=time.process_time()
print((et-st)*1e3)
print(data)
#print(y2.sum(),suma)

"""
cpixx.spline2.restype = np.ctypeslib.ndpointer(dtype=float)
cpixx.spline2.argtypes = [
    np.ctypeslib.ndpointer(dtype=float),
    np.ctypeslib.ndpointer(dtype=float),
    ctypes.c_int,
]
y2=cpixx.spline2(x,y,n)
print(y2)
#print(y2.sum())
"""
exit(0)

#########################################
#Test
#########################################
cpixx.test_cpixx()
exit(0)
