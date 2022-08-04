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
        Verbose.VERBOSITY=1
        print("Basic definition:")
        sg=Spangler(nspangles=3,body_hash="123",spangle_type=GRANULAR_SPANGLE)
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Equ->Obs:\n",sg.M_obs2ecl)
        print("Equ->Obs:\n",sg.M_equ2obs)

        print("\nAnother definition:")
        sg=Spangler(nspangles=3,body_hash="123",n_equ=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)
        print("Equ->Obs:\n",sg.M_equ2obs)

        print("\nDefinition observer:")
        sg=Spangler(nspangles=3,body_hash="123",n_equ=[0,1,0],n_obs=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)
        print("Equ->Obs:\n",sg.M_equ2obs)
        Verbose.VERBOSITY=0

    def test_pop(self):
        sg=Spangler(nspangles=100,body_hash="123")
        sg.populate_spangler(geometry="ring",scale=1,seed=1,gaps=[[0,0.2],[0.5,0.1]],boundary=0)
        print_df(sg.data.head(5))
        sg.update_positions(n_obs=[1,1,1])
        print_df(sg.data.head(5))
        sg=Spangler(nspangles=1000,body_hash="123",n_equ=[1,0,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        print(sg.nspangles,sg.sample.N,len(sg.data))
        
    def test_plot3d(self):
        Verbose.VERBOSITY=0
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        sg.plot3d(factor=1.3,c='b',s=3)
        
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",scale=2,seed=1,boundary=0)
        sg.plot3d(factor=0.4)
        #sg.populate_spangler(geometry="ring",scale=2,seed=1,gaps=[[0,0.2],[0.5,0.1]],boundary=0)
        Verbose.VERBOSITY=0

        
    def test_plotobs(self):
        Verbose.VERBOSITY=0
        sg=Spangler(nspangles=1000,body_hash="123",n_equ=[1,1,1])
        
        sg.populate_spangler(geometry="circle",scale=2,seed=1,boundary=0)
        sg.plot_obs()

        sg.populate_spangler(geometry="ring",scale=2,seed=1,gaps=[[0,0.2],[0.5,0.1]],boundary=0)
        sg.plot_obs()
        
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        sg.update_positions(n_obs=[1,0,1])
        sg.plot_obs()

        Verbose.VERBOSITY=0

        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
