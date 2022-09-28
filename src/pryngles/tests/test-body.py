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
        
        B=Body("Body",BODY_DEFAULTS,None,bhash='B',m=2,x=2)
        
        print(B)
        print(B.m)
        
        B.update_body(bhash="B")
        print(B)
        
        C=Body("Body",BODY_DEFAULTS,B,bhash="C")
        print(C)
        print(B)
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_spangle(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        #Create body
        B=Body("Body",BODY_DEFAULTS,None,bhash='B',m=2,x=2)
        B.spangle_body()
        B.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_star(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        S=Star()
        print(S)

        #Check derived properties
        self.assertEqual(np.isclose([S.wrot],
                                    [2*np.pi/BODY_DEFAULTS["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        
        S.update_star(m=2,limb_coeffs=[1,1])
        print(S)
        
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Star(primary="Nada"))     
        
        S=Star(nspangles=270,i=45*Consts.deg)
        S.spangle_body()
        
        print_df(S.sg.data.tail())
        
        S.sg.set_observer()
        S.sg.set_luz()
        S.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE

    def test_planet(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        S=Star()

        #Check exception: primary is mandatory for planets
        self.assertRaises(ValueError,lambda:Planet())

        P=Planet(primary=S)
        print(P.bhash)
        
        #Check derived properties
        self.assertEqual(np.isclose([P.wrot],
                                    [2*np.pi/BODY_DEFAULTS["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        
        #Check a non-existing property
        P.update_planet(vz=0.2)
        print(P)
        
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Planet(primary="Nada"))
        
        P.update_body(nspangles=250)
        P.spangle_body()
        print_df(P.sg.data.tail())
        
        P.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_ring(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        #Define first star and planet
        S=Star()
        P=Planet(primary=S)

        self.assertRaises(ValueError,lambda:Ring())
        R=Ring(primary=P)
        
        R.update_ring(fe=3)
        print(R)
        
        R.update_body(nspangles=250,i=60*Consts.deg,roll=0*Consts.deg)
        R.spangle_body()
        print_df(R.sg.data.tail())
        R.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
