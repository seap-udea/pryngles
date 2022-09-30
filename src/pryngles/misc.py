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
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: Miscellaneous

from pryngles import *

# ## External modules

from collections import OrderedDict as odict
from collections.abc import Iterable
import inspect

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
        """Print DataFrame.
        
        Parameters:
            df: Pandas DataFrame:
                DataFrame to print.
        """
        display(HTML(df.to_html()))
        
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
                
    def get_methods(my_class):
        """Get a list of the methods for class my_class
        """
        return sorted([member[0] for member in inspect.getmembers(my_class) if '__' not in member[0]])
                
Misc.__doc__=Misc_doc


