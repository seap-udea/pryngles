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
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,a=2,radius=2)
        M=sys.add("Planet",primary=P,a=2,radius=2)
        R=sys.add("Ring",primary=P,fi=1.3,fe=2.3)
        
        for particle in sys.sim.particles:
            print(particle)
            
        print(sys)
        print(sys.nbodies,sys.nparticles,sys.nsources)

        Verbose.VERBOSITY=VERB_NONE
                
    def test_system_remove(self):
        
        Verbose.VERBOSITY=VERB_ALL
        
        sys=System()
        S=sys.add(bhash="Star",m=8,radius=4,x=5,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",a=2,radius=2)
        M=sys.add("Planet",primary=P,bhash="Moon",a=2,radius=2)
        R=sys.add("Ring",primary=P,bhash="Ring",fi=1.3,fe=2.3)
        print(sys.bodies)
        sys.remove("Ring")
        print(sys.bodies)
        sys.remove("Planet")
        print(sys.bodies)
        sys.remove("Star")
        print(sys.bodies)

        Verbose.VERBOSITY=VERB_NONE

    def test_integrate(self):
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=False)
        S=sys.add(bhash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,a=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,a=1,M=120*Consts.deg)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        sys.sim.integrate(1)
        rb.OrbitPlot(sys.sim)
        
        Verbose.VERBOSITY=VERB_NONE

    def test_spangle(self):
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=False)
        S2=sys.add(bhash="Star2",nspangles=nspangles,m=8,radius=1,x=10,vy=-2)
        S1=sys.add(bhash="Star1",nspangles=nspangles,m=9,radius=1,x=0,vy=+2)
        P=sys.add("Planet",primary=S1,bhash="Planet",nspangles=nspangles,radius=0.2,a=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,a=1,M=120*Consts.deg)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        
        sys.spangle_system()
        
        #Check addition columns
        print(sys.nsources)
        print(sys.sources)
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
        S=sys.add(hash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,a=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,a=1,M=120*Consts.deg)
        R=sys.add("Ring",primary=P,bhash="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
        print(P.radius)
        sys.update_body(P,radius=0.5)
        print(P.radius)
        sys.update_body("Ring",fe=3.0)
        print(R.radius)
        sys.spangle_system()
        self.assertRaises(AssertionError,lambda:sys.update_body("Ring",fe=3.0))

        Verbose.VERBOSITY=VERB_NONE

    def test_obs(self):
        
        global sys
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        
        #Define system
        sys=System(resetable=True)
        
        #Add objects
        S=sys.add(bhash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,a=2)
        
        #Test setting observer without spangling
        self.assertRaises(AssertionError,lambda:sys.set_observer(nvec=[1,0,0]))
        
        #Spangle system
        sys.spangle_system()
        
        sys.sg.plot3d(center_at="Planet",not_plot=["Star"])
        
        sys.set_observer(nvec=[-1,0,0])
        sys.sg.plot3d()
        
        sys.set_observer(nvec=[0,0,1])
        sys.sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE

    def test_reset(self):
        
        Verbose.VERBOSITY=VERB_NONE
        
        nspangles=100
        sys=System(resetable=True)
        S=sys.add(bhash="Star",nspangles=nspangles,m=8,radius=1,x=0,vy=2)
        P=sys.add("Planet",primary=S,bhash="Planet",nspangles=nspangles,radius=0.2,a=2)
        M=sys.add("Planet",primary=P,bhash="Moon",nspangles=nspangles,radius=0.1,a=1,M=120*Consts.deg)
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


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
