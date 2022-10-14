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
    def test_Orbit(self):
        
        #Quintuple system
        S1=Orbit(m1=1,m2=1,a=1,e=0.7,M=0)
        S2=Orbit(m1=1,m2=1,a=1,e=0,M=0)
        S3=Orbit(S1,S2,a=5,e=0)
        S4=Orbit(S3,m2=1,a=20,e=0,E=45*Consts.deg)
        S4.calculate_orbit()
        print(S4.get_states())
        print(S4.Ps)
        
        #Using custom units
        #Quintuple system

        #Initialize positions
        hn=Orbit(
            m1=1,
            m2=Orbit(m1=1e-3,m2=1e-7,a=0.5,e=0.0,units=units),
            units=units,
            a=20,e=0.0)
        hn.calculate_orbit()
        sim,states=hn.get_states()
        print(states)
        
        #SImple 
        
        #Use this code to animate:
        #Plot.animate_rebound(S4.sim)
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
