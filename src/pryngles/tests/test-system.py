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
    def test_system_init(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        print(sys.nbodies)
        print(sys.sim.G)
        print(sys.ul,sys.um,sys.ut)
        
        sys=System(units=['m','kg','s'])
        print(sys.nbodies)
        print(sys.sim.G)
        print(sys.ul,sys.um,sys.ut)
        
        print(sys)

        Verbose.VERBOSITY=VERB_NONE

    def test_system_add(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,a=2,radius=2)
        M=sys.add("Planet",primary=P,a=2,radius=2)
        R=sys.add("Ring",primary=P,fi=1.3,fe=2.3)
        for particle in sys.sim.particles:
            print(particle)
            
        print(sys)

        Verbose.VERBOSITY=VERB_NONE
                
    def test_system_remove(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(hash="Star",m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,hash="Planet",a=2,radius=2)
        M=sys.add("Planet",primary=P,hash="Moon",a=2,radius=2)
        R=sys.add("Ring",primary=P,hash="Ring",fi=1.3,fe=2.3)
        print(sys.bodies)
        sys.remove("Ring")
        print(sys.bodies)
        sys.remove("Planet")
        print(sys.bodies)
        sys.remove("Star")
        print(sys.bodies)

        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
