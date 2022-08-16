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
    def test_init(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        #Define first star and planet
        S=Star()
        P=Planet(primary=S)

        self.assertRaises(ValueError,lambda:Ring())
        R=Ring(primary=P)
        
        print(R.physics)
        print(R.optics)
        print(R.hash)
        
        R.update_body(physics=dict(fe=3))
        print(R.physics)
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_sp(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        S=Star()
        P=Planet(primary=S)
        R=Ring(primary=P,physics=dict(fi=1.5,fe=5.0,i=30*Consts.deg,roll=45*Consts.deg))
        
        R.spangle_body(preset=True)
        R.sp.plot3d(factor=1.1)
        
        Verbose.VERBOSITY=VERB_NONE
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
