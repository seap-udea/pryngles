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

# # Pryngles module: physics
# 
# Template of a module

from pryngles import *

# ## External modules

#@external
import numpy as np
#@end

# ## Scatterer Class

# ### Docstring

Scatterer_doc="""This is the basic class of a scatterer
"""

# ### Class structure

class Scatterer(PrynglesCommon):
    def __init__(self):
        pass
    
Scatterer.__doc__=Scatterer_doc

#Define system
nspangles=1000
sys=System()
S=sys.add(nspangles=nspangles)
P=sys.add("Planet",primary=S,nspangles=nspangles,radius=0.2,x=2)
sys.spangle_system()

sys.sg.plot2d(include=["Planet"])

sys.update_perspective(n_obs=sci.direction(90,0))

sys.sg.plot3d()

sys.sg.plot2d()

get_ipython().run_line_magic('pinfo', 'sys.sg.plot2d')

sys.sg.plot2d(show_azim=True)

#@test:template
#@end

# ### Method development

#@method:Foo
def method2(self):
    return self.a+self.b
#@end

Foo.method2=method2

#@test:template
#@end

