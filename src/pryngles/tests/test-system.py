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
        
        sys=System(resetable=True)
        print("Nbodies = ",sys.nbodies)
        print("G constant = ",sys.sim.G)
        print("G constant = ",sys.units)
        print("Canonical units = ",sys.ul,sys.um,sys.ut)
        
        sys=System(units=['m','kg','s'])
        print("Nbodies = ",sys.nbodies)
        print("G constant = ",sys.sim.G)
        print("G constant = ",sys.units)
        print("Canonical units = ",sys.ul,sys.um,sys.ut)
        print(sys)
        
        sys.save_to("/tmp/system.pkl")
        print(sys.sim.status())
        sys2=System("/tmp/system.pkl")
        print(sys2.sim.status())
        
        print(sys.sim.N)

        Verbose.VERBOSITY=VERB_NONE

    def test_system_add(self):
        global sys
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,radius=2,x=1)
        M=sys.add("Planet",bhash="Moon",primary=P,radius=2,x=2)
        R=sys.add("Ring",primary=P,fi=1.3,fe=2.3)
        
        for particle in sys.sim.particles:
            print(particle)
            
        print(sys)
        print(sys.nbodies,sys.nparticles)

        Verbose.VERBOSITY=VERB_NONE
                
    def test_system_remove(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(bhash="Star",m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",radius=2,x=2)
        M=sys.add("Planet",primary=P,bhash="Moon",radius=2,x=2)
        R=sys.add("Ring",primary=P,bhash="Ring",fi=1.3,fe=2.3)
        print(sys.bodies)
        sys.remove("Ring")
        print(sys.bodies)
        sys.remove("Planet")
        print(sys.bodies)
        sys.remove("Star")
        print(sys.bodies)

        Verbose.VERBOSITY=VERB_NONE

    def test_obs(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        
        #Define system
        sys=System(resetable=True)
        
        #Add objects
        S=sys.add(nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,nspangles=nspangles,radius=0.2,x=2)
        
        #Test setting observer without spangling
        self.assertRaises(AssertionError,lambda:sys.set_observer(nvec=[1,0,0]))
        
        #Spangle system
        sys.spangle_system()
        
        sys.set_observer(nvec=[-1,0,0])
        sys.sg.plot3d()
        
        sys.set_observer(nvec=[0,0,1])
        sys.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE

    def test_setluz(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        nspangles=500
        sys=System()
        S=sys.add("Star",nspangles=nspangles,m=1,radius=1,x=0,y=0,vy=0)
        D=sys.add("Ring",bhash="Disk",primary=S,nspangles=nspangles,fi=5,fe=10,i=0*Consts.deg)
        P=sys.add("Planet",primary=S,nspangles=nspangles,radius=0.2,m=1e-3,x=2)
        R=sys.add("Ring",primary=P,nspangles=nspangles,fi=1.5,fe=2.0,i=-20*Consts.deg)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,m=1e-6,x=2.5,y=-0.2)
        K=sys.add("Ring",bhash="Cronoring",primary=M,nspangles=nspangles,fi=1.1,fe=1.5,i=20*Consts.deg)
        sys.spangle_system()
        sys.sg.plot3d(center_at="Ring",not_plot=["Disk"])
        cond=(sys.sg.data.bhash=="Moon")&(sys.sg.data.hidden_by_luz!="")
        print_df(sys.sg.data.loc[cond,["hidden_by_luz"]].head(10))

        Verbose.VERBOSITY=VERB_NONE

    def test_spangle(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=False)
        S2=sys.add(bhash="Star2",nspangles=nspangles,m=8,radius=1,x=10,vy=-2)
        S1=sys.add(bhash="Star1",nspangles=nspangles,m=9,radius=1,x=0,vy=+2)
        P=sys.add("Planet",primary=S1,bhash="Planet",nspangles=nspangles,radius=0.2,x=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,x=3)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        
        sys.spangle_system()
        
        #Check addition columns
        print(sys.source)
        print(sys.sg.data.columns)
        
        #Check save
        sys.save_to("/tmp/system.pkl")
        
        #Check plot
        #sys.sp.plot3d(center_at="Ring",not_plot=["Star1","Star2"])
        sys.sg.plot3d()

        Verbose.VERBOSITY=VERB_NONE

    def test_update(self):
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=True)
        S=sys.add("Star",bhash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,x=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,x=3)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        print(P.radius)
        sys.update_body(P,radius=0.5)
        print(P.radius)
        sys.update_body("Ring",fe=3.0)
        print(R.radius)
        sys.spangle_system()
        self.assertRaises(AssertionError,lambda:sys.update_body("Ring",fe=3.0))

        Verbose.VERBOSITY=VERB_NONE

    def test_reset(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=True)
        S=sys.add("Star",bhash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,x=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,x=3)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        sys.spangle_system()

        #All transformations from here are not stored
        sys.sg.plot3d()
        sys.set_observer(nvec=[0,0,-1])
        sys.sg.plot3d()

        #All transformations from here are not stored
        sys.reset()
        sys.sg.plot3d()

        Verbose.VERBOSITY=VERB_NONE

    def test_int(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System()
        S=sys.add("Star",bhash="Star",nspangles=nspangles,m=1,radius=1,x=0,y=0,vy=0)
        M=sys.add("Planet",primary=S,bhash="Moon",nspangles=nspangles,radius=0.1,m=1e-6,
                  x=0,y=-3,vx=np.sqrt(sys.sim.G/3),vy=0)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,m=1e-3,
                  x=2,vy=np.sqrt(sys.sim.G/2))
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=20*Consts.deg)
        sys.spangle_system()

        sys.integrate(1)
        
        sys.set_observer([0,0,1])
        sys.set_luz()
        
        sys.integrate(1)
        
        sys.set_observer([0,0,1])
        sys.set_luz()

        sys.sg.plot3d()
        sys.sg.plot3d(center_at="Ring")
        
        Verbose.VERBOSITY=VERB_NONE

    def test_anim(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System()
        S=sys.add("Star",bhash="Star",nspangles=nspangles,m=1,radius=1,x=0,y=0,vy=0)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,m=1e-3,
                  x=2,vy=np.sqrt(sys.sim.G/2))
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,m=1e-6,
                  x=0,y=-3,vx=np.sqrt(sys.sim.G/3),vy=0)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=20*Consts.deg)
        sys.spangle_system()

        sys.animate_integration(filename="/tmp/pryngles-dynamics.gif",tini=0,tend=10,interval=1000,nsnap=5,coords="ecl")
        
        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
