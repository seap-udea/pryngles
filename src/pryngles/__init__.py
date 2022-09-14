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
##################################################
# PRELIMINARY INIT COMMANDS
##################################################
#--END OF TEMPLATE--#

#!/usr/bin/env python
# coding: utf-8

# # PlanetaRY spanGLES

# This is the initialization file of the `Pryngles` package.

# ## Warnings

import unittest
import warnings
import dill
import inspect
from collections import OrderedDict as odict
from collections.abc import Iterable
from copy import deepcopy
from colorsys import hls_to_rgb
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ## Jupyter compatibilty

"""
The purpose of the get_ipython class is to provide some response in the python 
script resulting from the conversion of this notebook.

If you want to add another IPyhton function resulting from a magic command to the class, 
please verify in the resulting python script the corresponding IPython command.

For instance, the magic "%matplotlib nbagg" is converted into:

    get_ipython().magic('matplotlib nbagg',globals())

So, the method "magic" should be add to the get_ipython() class.        
"""
from IPython.display import HTML, Image, display
import IPython.core.autocall as autocall
from IPython import get_ipython
import sys

try:
    cfg=get_ipython().config
except AttributeError:
    def Image(url="",filename="",f=""):
        pass
    class get_ipython(object):
        def run_line_magic(self,*args):
            pass
        def run_cell_magic(self,x,y,z):
            pass
        def magic(self,command,scope=globals()):
            pass

#Magics can only be located from here
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

# ## Verbosity

"""Verbosity levels:
SIMPLE: Simple messages.
SYSTEM: System operations.
VERIFY: Message to verify operations
"""
VERB_NONE=0
VERB_SIMPLE=1
VERB_SYSTEM=2
VERB_VERIFY=3
VERB_DEEP=4
VERB_ALL=100

class Verbose(object):
    """Verbose print in the package
    
    Attributes:
        VERBOSITY: int, default = 0:
            Level of verbosity.
            
    Methods:
        print(level,msg):
            Print a message if level<=VERBOSITY.
    
    Example:
    
        Verbose.print(1,"Hello world") #No print
        
        Verbose.print(0,"Hello world") #Print

        Verbose.VERBOSITY=1
        Verbose.print(1,"Hello world") #Print
        
        Verbose.VERBOSITY=2
        Verbose.print(1,"Hello world") #Print
        
        Verbose.VERBOSITY=2
        Verbose.print(4,"Hello world") #No print
    """
    VERBOSITY=VERB_ALL
    def print(level,*args):
        if level<=Verbose.VERBOSITY:
            print("  "*level+f"VERB{level}::{inspect.stack()[1][3]}::",*args)

verbose=Verbose.print

# ## PrynglesCommon
# 
# Many of the classes in Pryngles inherite methods of this common class

class PrynglesCommon(object):
    
    def save_to(self,filename):
        """Save object to a binary file
        
        Parameters:
            filename: string:
                Name of the file where the object will be stored.
                
            compressed: boolean, default = False:
                If True the file will be stored compressed.
        
        Notes:
            Based on https://betterprogramming.pub/load-fast-load-big-with-compressed-pickles-5f311584507e.
        """
        verbose(VERB_SYSTEM,f"Saving object to {filename}")
        pikd = open(filename,"wb")
        dill.dump(self, pikd)
        pikd.close()
            
    def load_from(self,filename,compressed=False):
        verbose(VERB_SYSTEM,f"Loading object from {filename}")
        pikd = open(filename,"rb")
        data = dill.load(pikd)
        pikd.close()
        verbose(VERB_VERIFY,f"Transferring data to new object")
        self.__dict__=data.__dict__
        return data
    
    def __str__(self):
        #Remove private attributes
        return str({k:v for k,v in self.__dict__.items() if k[0]!='_'})

# ## Miscelaneous Class

Misc_doc="""
Miscelaneous routines.

This is a set of util routines intended for a diversity of purposes.

Routines included:

    get_data(file)
""";

class Misc(object):
    def get_data(path):
        """
        Get the full path of the `datafile` which is one of the datafiles provided with the package.
        
        Parameters:
            datafile: Name of the data file, string.
            
        Return:
            Full path to package datafile in the python environment.
            
        """
        return os.path.join(ROOTDIR,'data',path);
    
    def print_df(df):
        """
        Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        display(HTML(df.to_html()))
        
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
                c=Circle((L,S),dL/2.5,color=Misc.rgb([H,L,S]))
                ax.add_patch(c)
                ax.text(L,S,f"S={S:.1g},L={L:.1g}",ha='center',va='center',fontsize=6,color='y')
        ax.axis("off")
        ax.axis("equal")
        plt.tight_layout()
            
Misc.__doc__=Misc_doc

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
Misc.flatten=flatten

# ## Pryngles modules

from pryngles.version import *
from pryngles.consts import *
from pryngles.science import *
from pryngles.plot import *
from pryngles.sampler import *
from pryngles.spangler import *
from pryngles.body import *
from pryngles.system import *
from pryngles.legacy import *
#Reset verbosity
Verbose.VERBOSITY=VERB_NONE

# ## Tests



