##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#

# PlanetaRY spanGLES                                             #
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
# Main contributors:                                             #
#   Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado         #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
	def test_scatterer(self):
	    Verbose.VERBOSITY=VERB_ALL
	    
	    au = 1.496e+11
	    sys=System(units=["au","msun","yr"])
	    S=sys.add(radius=Consts.rsun/au)
	    P=sys.add("Planet",primary=S,radius=Consts.rsaturn/au,a=3,nspangles=1000)
	    sys.initialize_simulation()
	    sys.spangle_system()
	
	    incli = 0
	    azim = 179
	    sys.update_perspective(n_obs=Science.direction(azim,incli))
	    
	    # Test Scatterer class
	    test = Scatterer(fname_planet="./data/fou_lambert.dat",fname_ring=None)
	    out_dict = test.scattering(sys,normalize=True)
	    F = test.lambertian_test(np.arccos(test.phase_angle))/np.pi
	    print("Flux from pryngles: ", out_dict["Stot"][0])
	    print("Flux from theory: ", F)
	
	    Verbose.VERBOSITY=VERB_NONE
	
if __name__=="__main__":
    unittest.main(argv=['first-arg-is-ignored'],exit=False)
    