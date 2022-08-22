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
    def test_temp(self):
        
        #Test
        Science.template()
        """
        self.assertEqual(self.P.Nr,8,True)
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        
    def test_coords(self):
        
        #Test spherical
        rqf=Science.spherical([1,1,0])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,+1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,-1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,-1,1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,+1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([-1,-1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)   
        rqf=Science.spherical([+1,-1,-1])
        print(rqf[0],rqf[1]*Consts.rad,rqf[2]*Consts.rad)
        
        #Test cartesian
        xyz=Science.cartesian([1,0,0])
        print(xyz) 
        xyz=Science.cartesian([1,45*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,135*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,225*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,315*Consts.deg,45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,45*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,135*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,225*Consts.deg,-45*Consts.deg])
        print(xyz) 
        xyz=Science.cartesian([1,315*Consts.deg,-45*Consts.deg])
        print(xyz) 

    def test_rot(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        #Test rotation
        Msys2uni,Muni2sys=Science.rotation_matrix([0,0,1],0)
        print(Msys2uni)

        Msys2uni,Muni2sys=Science.rotation_matrix([0,0,-1],0)
        print(Msys2uni)

        Verbose.VERBOSITY=VERB_NONE

    def test_limb(self):
        
        Verbose.VERBOSITY=VERB_ALL

        cs=[np.random.rand()]
        I=Science.limb_darkening(0.8,cs)
        print(I)
        
        fig=plt.figure()
        ax=fig.gca()

        rhos=np.linspace(0,1,100)
        coefs=[0.6550]
        ax.plot(rhos,Science.limb_darkening(rhos,coefs))
        coefs=[0.6022,0.0654]
        ax.plot(rhos,Science.limb_darkening(rhos,coefs))
        coefs=[0.9724,-0.4962,0.2029]
        ax.plot(rhos,Science.limb_darkening(rhos,coefs))    
        coefs=[-0.2018,2.1000,-2.0247,0.7567]
        ax.plot(rhos,Science.limb_darkening(rhos,coefs))
        Plot.pryngles_mark(ax)
        
        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
