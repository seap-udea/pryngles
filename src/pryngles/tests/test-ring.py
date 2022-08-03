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
        R=Ring(primary=P,optics=dict(albedo_gray_normal=0.5,tau_gray_optical=0.2,nspangles=100))  
        R.spangle_body()
        print(R.sp.N)
        print(R.optics)
        
        Misc.print_html(R.sg.df.head(5).to_html())
        
        spangled=None
        spangled=dict(color='r')
        R.sp.plot(spangled=spangled)
        #"""
        #print(len(R.spangles))
        print(R.sp.aes)
        
    def test_plot(self):
        S=Star()
        P=Planet(primary=S)
        R=Ring(primary=P,optics=dict(nspangles=100))  
        R.spangle_body()
        R.plot_body()
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
