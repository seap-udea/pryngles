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
import unittest
from pryngles import *
class Test(unittest.TestCase):
    def test_fun(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        B=Body(BodyDefaults,"Body",None,dict(x=0),dict(mass=2.3),dict())
        
        print(B)
        print(B.orbit)
        print(B.physics)
        print(B.optics)
        
        B.update_body(orbit=dict(m=2))
        print(B.orbit)
        
        C=Body(BodyDefaults,"Body",B,dict(),dict(),dict())
        print(C)
        print(B)
        
        Verbose.VERBOSITY=VERB_NONE
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
