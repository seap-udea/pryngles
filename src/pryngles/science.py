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

# # Pryngles module: Science

from pryngles import *

# ## External modules

import numpy as np

# ## The Science class
# 
# The Science class is a class with routines intended to perform a wide diversity of mathematical, physical and astronomical calculations.

class Science(object):
    pass

def xyz2rqf(xyz):
    """
    Transform cartesian coordinates into spherical coordinates

    Parameters:

        xyz: array (3):
            Cartesian coordinates

    Return:

        rqf: array (3):
            Spherical coordinates (r, theta, phi) where theta is azimutal angle and phi is 
            elevation (complement of polar angle).                

            Notice that this convention is different than that of regular vectorial calculus
            where spherical coordinates are (r,theta,phi), but theta is the polar angle and phi 
            the ezimutal one.

    """
    rho=(xyz[0]**2+xyz[1]**2)**0.5
    r=(rho**2+xyz[2]**2)**0.5
    theta=np.arctan2(xyz[1],xyz[0])
    phi=np.arctan2(xyz[2],rho)
    return np.array([r,theta,phi])

Science.xyz2rqf=xyz2rqf


