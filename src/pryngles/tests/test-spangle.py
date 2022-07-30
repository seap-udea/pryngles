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
    def test_sample(self):

        #Generate rings
        sp=Spangle()
        print("Default : ",sp)
        
        #Set positions
        sp.set_position(
            [[1,2,3],[4,5,6]],
            [[7,8,9],[10,11,12]]
        )
        
        sp.set_orientation([[1,1,1],[1,1,1]])
        
        #Set albedo
        sp.set_optical(albedo_gray_normal=0.5)

        sp.set_optical(tau_gray_optical=5)

        print("Modified : ",sp)

        #Errors
        self.assertRaises(KeyError,lambda:sp.set_optical(albedo=2))
        """
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
