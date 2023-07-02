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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External required packages
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import unittest
import warnings
import dill
import inspect
import sigfig

from copy import deepcopy
from collections import OrderedDict as odict

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constants and Configurations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Version
from pryngles.version import *

# Ignore warnings
warnings.filterwarnings('ignore')

# Verbosity levels: see help(Verbose)
VERB_NONE=0
VERB_SIMPLE=1
VERB_SYSTEM=2
VERB_VERIFY=3
VERB_DEEP=4
VERB_ALL=100

# Type of figure output: 'save', 'plot', None
VERB_TEST_FIG = 'save'
VERB_TEST_DIR = '/tmp'

# Pryngles Data Repository
DATA_INDEX=dict(
    baseurl = 'https://docs.google.com/feeds/download/spreadsheets/Export?exportFormat=csv&key=',
    downurl = 'https://docs.google.com/uc?export=download&id=',
    filename = 'pryngles_data_index.csv',
    fileid = '17rgmzuENn_5jEO2rzDzngJLpZT1P9rPsQr3PNQGFkGs',
    col_filename = 'filename',
    col_fileid = 'fileid',
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Common classes of all package
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Verbose(object):
    """Verbose print in the package
    
    Attributes:
        VERBOSITY: int, default = 0:
            Level of verbosity.
            
            Verbosity levels:
                SIMPLE: Simple messages.
                SYSTEM: System operations.
                VERIFY: Message to verify operations
                DEEP: Deep debugging messages
                ALL: All debugging messages
                
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

    def save_test_fig(suffix):
        if VERB_TEST_FIG == 'save':
            for ifig in plt.get_fignums():
                fig = plt.figure(ifig)
                fig.savefig(f'{VERB_TEST_DIR}/test-{suffix}-{ifig:03d}.png')
        elif VERB_TEST_FIG == 'show':
            plt.show()
        else:
            pass
        
"""Mandatory alias

The method `Verbose.print` will be used across many modules in
pryngles, so it shouldn't be removed from here
"""
verbose=Verbose.print
            
class PrynglesCommon(object):
    """Base class of the package.
    
    All major classes are children of PrynglesCommon class.
    """
    def __init__(self):
        pass
    
    def save_to(self,filename):
        """Save object to a binary file
        
        Parameters:
            filename: string:
                Name of the file where the object will be stored.
        
        Notes:
            Based on https://betterprogramming.pub/load-fast-load-big-with-compressed-pickles-5f311584507e.
        """
        verbose(VERB_SYSTEM,f"Saving object to {filename}")
        pikd = open(filename,"wb")
        dill.dump(self, pikd)
        pikd.close()
            
    def load_from(self,filename):
        """Read object from a binary file.
        
        Parameters:
            filename: string:
                Name of the file where the object is stored.        
        """
        verbose(VERB_SYSTEM,f"Loading object from {filename}")
        pikd = open(filename,"rb")
        data = dill.load(pikd)
        pikd.close()
        verbose(VERB_VERIFY,f"Transferring data to new object")
        self.__dict__=data.__dict__
        return data
    
    def __str__(self):
        """Show content of an object
        
        This method determines the default behavior of the command:
        
            print(object)
        """
        #Remove private attributes
        return str({k:v for k,v in self.__dict__.items() if k[0]!='_'})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load all modules of the package
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spin-off
from pryngles.hinb import * 
from pryngles.fibosampler import * 

#Utility modules 
from pryngles.util import *

#Physical modules 
from pryngles.optics import * 

#Legacy module 
from pryngles.legacy import *

#Core modules
from pryngles.spangler import * 
from pryngles.system import *

#Photometry
from pryngles.photometry import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Aliases and initialization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Reset verbosity
Verbose.VERBOSITY=VERB_NONE

#This aliases does not work in modules
print_df=Misc.print_df
sci=Science

