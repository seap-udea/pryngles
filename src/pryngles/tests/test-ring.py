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
        
        #Define first star and planet
        S=Star()
        P=Planet(primary=S)

        self.assertRaises(ValueError,lambda:Ring())
        R=Ring(primary=P)
        
        print(R.physics)
        print(R.optics)
        print(R.hash)
        
    def test_update(self):

        #Define first star and planet
        S=Star()
        P=Planet(primary=S)
        R=Ring(primary=P)
        
        R.update_body(physics=dict(fe=3))
        print(R.physics)
        
        #Check derived properties
        print(R.M_equ2ecl)
        print(R.M_ecl2equ)
        print(R.nr_ecl)
        print(R.nr_equ)
        
    def test_spangle(self):
        S=Star()
        P=Planet(primary=S)
        R=Ring(primary=P)  
        R.spangle_body()
        print(R.sp.N)
        print(R.optics)
        print(R.spangles[0].xyz)
        
        #return
        #"""
        spangled=None
        spangled=dict(color='r')
        R.sp.plot(spangled=spangled)
        #"""
        #print(len(R.spangles))
        print(R.sp.aes)
        print(R.spangles[0].asp)
        
    def test_plot(self):
        S=Star()
        P=Planet(primary=S)
        R=Ring(primary=P)  
        R.spangle_body()
        print(R.spangles[0].xyz)
        #R.plot_body()
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
