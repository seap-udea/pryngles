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
        f=Foo()
        print(f)
        
    def test_fun(self):
        f=Foo(-1,5)
        c=f.method2()
        print(c)
        
        """
        #Type of tests
        self.assertRaises(ValueError,lambda:Planet())
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        """
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
