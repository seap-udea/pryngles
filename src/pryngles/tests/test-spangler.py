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
        
        Verbose.VERBOSITY=VERB_ALL
        
        print("Basic definition:")
        sg=Spangler(nspangles=3,spangle_type=GRANULAR_SPANGLE)
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Equ->Obs:\n",sg.M_obs2ecl)

        print("\nAnother definition:")
        sg=Spangler(nspangles=3,n_equ=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)

        print("\nDefinition observer:")
        sg=Spangler(nspangles=3,body_hash="123",n_equ=[0,1,0],n_obs=[1,0,0])
        Misc.print_df(sg.data.head(5))
        print("Equ->Ecl:\n",sg.M_equ2ecl)
        print("Obs->Ecl:\n",sg.M_obs2ecl)
        
        Verbose.VERBOSITY=VERB_NONE

    def test_pop(self):
        Verbose.VERBOSITY=VERB_ALL
        
        #Using preset
        sg=Spangler(nspangles=850)
        sg.populate_spangler(geometry="ring",preset=True,ri=0.5)
        sg.sample.plot()
        print(sg.sample.N)
        print(sg.nspangles)
    
        #Ring
        sg=Spangler(nspangles=100)
        sg.populate_spangler(geometry="ring",scale=1,seed=1,boundary=0)
        print_df(sg.data.head(5))
        sg.set_observer(nvec=[1,1,1])
        print_df(sg.data.head(5))
        
        #Sphere
        sg=Spangler(nspangles=1000,body_hash="123",n_equ=[1,0,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        print(sg.nspangles,sg.sample.N,len(sg.data))
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_plot3d(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        plt.close("all")
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",scale=2,seed=1)
        sg.set_luz(nvec=[1,1,0])
        sg.plot3d(factor=1.5)
        
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2,seed=1,boundary=0)
        sg.set_luz(nvec=[-1,1,-1])
        sg.set_observer(nvec=[1,0,1])
        sg.plot3d(show_hidden=True,factor=0.5,center_at="1234")

        Verbose.VERBOSITY=VERB_NONE

    def test_plotobs(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sg=Spangler(nspangles=2500,body_hash="123",n_equ=[1,1,1],center_ecl=[0,0,2])
        sg.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=2,seed=1,preset=True)
        sg.set_observer(nvec=[1,0,0])
        sg.set_luz(nvec=[1,1,1])
        sg.plot_obs()
        
        sg=Spangler(nspangles=500,body_hash="123",n_equ=[1,1,1],center_ecl=[1,1,1])

        sg.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=2,seed=1,preset=True)
        sg.set_observer(nvec=[1,0,0])
        sg.plot_obs()

        sg.populate_spangler(geometry="circle",spangle_type=GRANULAR_SPANGLE,scale=2,seed=1,boundary=0)
        sg.plot_obs()

        sg.set_luz(nvec=[0,0,-1])
        sg.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2,seed=1,ri=0.2,boundary=0)
        sg.plot_obs(center_at="123",show_hidden=0)
        
        Verbose.VERBOSITY=VERB_NONE

    def test_join(self):
        
        Verbose.VERBOSITY=VERB_ALL

        sg1=Spangler(nspangles=1000,body_hash="Ring",n_equ=[1,0,5])
        sg1.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)

        sg2=Spangler(nspangles=1000,body_hash="Planet",n_equ=[0,0,1])
        sg2.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=1,seed=1,preset=True)

        sgj=Spangler(spanglers=[sg1,sg2],n_obs=[1,0,0],n_luz=[+1]*3)

        sgj.plot3d(not_plot=["Ring"])
        sgj.plot_obs(show_hidden=0,not_plot=["Ring"])
        sgj.plot_obs(show_hidden=0,not_plot=[])
        
        Verbose.VERBOSITY=VERB_NONE

    def test_scale(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(center_ecl=[1,1,1],center_equ=[1,1,1])
        print_df(sg.data)

        sg.set_scale(5)
        print_df(sg.data)

        Verbose.VERBOSITY=VERB_NONE

    def test_reset(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(nspangles=100)
        sg.reset_state()
        print_df(sg.data.head())

        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
