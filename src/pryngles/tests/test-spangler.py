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
        sg=Spangler(nspangles=1,center_equ=[0,0,0],n_equ=[1,0,0])
        Misc.print_df(sg.data.head(1))

        print("\nCenter equ:")
        sg=Spangler(nspangles=3,center_equ=[0,0,1],n_equ=[0,1,0])
        Misc.print_df(sg.data.head(1))

        print("\nCenter ecl:")
        sg=Spangler(nspangles=3,center_ecl=[0,0,1],n_equ=[0,0,1])
        Misc.print_df(sg.data.head(1))

        print("\nRotation:")
        sg=Spangler(nspangles=3,w=30*Consts.deg,q0=40*Consts.deg,n_equ=[0,1,1])
        sg.set_positions(t=1)
        Misc.print_df(sg.data.head(1))

        print("\nJoin:")
        sg1=Spangler(sphash="Body 1",nspangles=3,w=40*Consts.deg,n_equ=[1,1,0])
        sg2=Spangler(sphash="Body 2",nspangles=3,w=30*Consts.deg,n_equ=[1,0,1])
        sg=Spangler(spanglers=[sg1,sg2])
        sg.set_positions(t=1)
        Misc.print_df(sg.data)

        Verbose.VERBOSITY=VERB_NONE

    def test_reset(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(nspangles=100)
        sg.reset_state()
        print_df(sg.data[["unset"]+list(SPANGLER_VISIBILITY_STATES)+list(SPANGLER_SOURCE_STATES)].head())

        Verbose.VERBOSITY=VERB_NONE

    def test_scale(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(center_ecl=[1,1,1],center_equ=[1,1,1])
        print_df(sg.data)

        sg.set_scale(5)
        print_df(sg.data)

        Verbose.VERBOSITY=VERB_NONE

    def test_pop(self):
        Verbose.VERBOSITY=VERB_ALL

        #No preset
        sg=Spangler(nspangles=850,n_equ=[1,0,0])
        sg.populate_spangler(shape="ring",
                             spangle_type=SPANGLE_GASEOUS,
                             scale=2,seed=1,ri=0.2)
        sg.sample.plot()
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        print_df(sg.data.head(3))

        #Using preset
        sg=Spangler(nspangles=850)
        sg.populate_spangler(shape="ring",
                             preset=True,
                             spangle_type=SPANGLE_SOLID_ROCK,ri=0.2)
        sg.sample.plot()
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        print_df(sg.data.head(3))
    
        #Sphere
        sg=Spangler(nspangles=100)
        sg.populate_spangler(shape="sphere",scale=3,seed=1,preset=True)
        sg.sample.plot(spangled=dict(color='r',alpha=0.1))
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        
        print_df(sg.data.head(3))
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_plot3d(self):
        global sg
        Verbose.VERBOSITY=VERB_SIMPLE

        #Sphere
        sg=Spangler(nspangles=100)
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ICE,preset=0,scale=3)
        sg.reset_state()
        
        cond=sg.data.z_ecl>0
        sg.data.loc[cond,"illuminated"]=True
        cond=sg.data.x_ecl>0
        sg.data.loc[cond,"visible"]=True
        cond=sg.data.y_ecl>0
        sg.data.loc[cond,"shadow"]=True
        cond=sg.data.f_equ>45*Consts.deg
        sg.data.loc[cond,"transmit"]=True

        sg.plot3d(statemark=0.5,coords="ecl")
    
        #No preset
        sg=Spangler(nspangles=850,n_equ=[1,1,1])
        sg.populate_spangler(shape="ring",preset=True,
                             spangle_type=SPANGLE_GRANULAR,
                             scale=2,ri=0.2)
        sg.data.illuminated=True
        sg.data.illuminated=True
        cond=sg.data.x_ecl>0
        sg.data.loc[cond,"visible"]=True
        sg.data.loc[cond,"transmit"]=True
        sg.plot3d(statemark=0.1)
    
        Verbose.VERBOSITY=VERB_NONE
        
    def test_setint(self):
        global sg
        
        Verbose.VERBOSITY=VERB_SIMPLE

        #No preset
        sg=Spangler(nspangles=50,sphash="Ring")
        sg.populate_spangler(shape="ring",seed=1,
                             spangle_type=SPANGLE_GRANULAR,
                             scale=2,ri=0.2)
        sg.data.illuminated=True
        sg.data.visible=True
        
        cond,n_int,d_int=sg.set_intersect(nvec=[1,0,1],center=[0,0,-1],
                                          sphash="Ring")
        sg._calc_qhulls()
        sg._plot_qhulls()

        #Plot 3d
        sg.plot3d(coords="int")
        plane=sg.qhulls["Ring"][0]["plane"]
        plane.plot_plane(ax=sg.ax3d,color='c',alpha=0.5)

        #Hulls
        print(sg.qhulls)
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_setobsluz(self):
        global sg
        
        Verbose.VERBOSITY=VERB_SIMPLE

        #Normal
        nspangles=10
        sg=Spangler(nspangles=nspangles,n_equ=[1,0,1],sphash="Planet")
        sg.populate_spangler(shape="sphere",preset=0,
                             spangle_type=SPANGLE_SOLID_ROCK,
                             scale=2)

        print_df(sg.data.loc[~sg.data.hidden,SPANGLER_KEY_FIELDS])

        sg.set_observer(nvec=[0,0,+1],center=None)
        sg.set_luz(nvec=[+1,0,0],center=None)

        sg.plot3d(coords="obs",statemark=1)

        #Semitransparent
        nspangles=50
        sg=Spangler(nspangles=nspangles,n_equ=[1,0,1],sphash="Planet")
        sg.populate_spangler(shape="sphere",preset=0,
                             spangle_type=SPANGLE_GASEOUS,
                             scale=2)
        sg.set_observer(nvec=[0,0,+1],center=None)
        sg.set_luz(nvec=[+1,0,0],center=None)
        sg.plot3d(statemark=1)

        Verbose.VERBOSITY=VERB_NONE
        
    def test_simplevis(self):
        global sg
        
        Verbose.VERBOSITY=VERB_SIMPLE

        plt.close("all")
        #Ring with semitransparent spangle: all illuminated, all visible, no transmission
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(shape="ring",ri=0.3,spangle_type=SPANGLE_GRANULAR,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([1,1,-1])
        sg.plot3d()

        #Ring with semitransparent spangle: all illuminated, all visible, no transmission
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(shape="ring",ri=0.3,spangle_type=SPANGLE_GRANULAR,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([-1,-1,-1])
        sg.plot3d()
        
        #Sphere with solid spangle: only illuminated 
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,preset=True,scale=3)
        sg.set_observer([1,0,1])
        sg.set_luz([0,0,1])
        sg.plot3d()
        
        #Sphere with stellar spangle: all illuminated, not all visible
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_STELLAR,preset=True,scale=3)
        sg.set_observer([1,0,1])
        sg.set_luz([0,0,1])
        sg.plot3d()

        #Sphere with semitransparent spangle: all illuminated, all visible
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_GASEOUS,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([1,0,0])
        sg.plot3d()

        #Two spheres
        sg1=Spangler(sphash="Planet 1",nspangles=100,center_equ=[-5,0,0])
        sg1.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ICE,preset=True,scale=3)
        
        sg2=Spangler(sphash="Planet 2",nspangles=100,center_equ=[+5,0,0])
        sg2.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,preset=True,scale=3)
        
        sg=Spangler(spanglers=[sg1,sg2])

        sg.set_observer([0,1,0])
        sg.set_luz(nvec=[1,0,0],center=[0,0,0],sphash="Planet 1")
        sg.set_luz(nvec=[-1,0,0],sphash="Planet 2")
        
        sg.plot3d()
        return
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_plot2d(self):
        
        global sg
        
        Verbose.VERBOSITY=VERB_SIMPLE
        
        plt.close("all")
        sg=Spangler(nspangles=2500,sphash="123",n_equ=[1,1,1],center_ecl=[0,0,2])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=2,seed=1,preset=True)
        #sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_GASEOUS,scale=2,seed=1,preset=True)

        sg.set_observer(nvec=[1,0,0])
        sg.set_luz(nvec=[1,1,1])
        fs=3
        sg.plot3d(coords="ecl")
        sg.plot2d(coords="ecl",fsize=fs)
        sg.plot2d(coords="luz",fsize=fs)
        sg.plot2d(coords="obs",fsize=fs)
        
        sg=Spangler(nspangles=200,sphash="123",n_equ=[1,1,1],center_ecl=[0,0,2])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=2,seed=1,preset=True)
        sg.set_observer(nvec=[1,0,0])
        sg.set_luz(nvec=[1,1,1])
        sg.plot2d(coords="obs",show_azim=True,fsize=5)
        
        Verbose.VERBOSITY=VERB_NONE

    def test_join(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        sg1=Spangler(nspangles=1000,sphash="Ring",n_equ=[1,0,5])
        sg1.populate_spangler(shape="ring",spangle_type=SPANGLE_GRANULAR,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)

        sg2=Spangler(nspangles=1000,sphash="Planet",n_equ=[0,0,1])
        sg2.populate_spangler(shape="sphere",spangle_type=SPANGLE_ATMOSPHERIC,scale=1,seed=1,preset=True)

        sgj=Spangler(spanglers=[sg1,sg2])
        
        sgj.set_observer([1,0,0.1])
        sgj.set_luz([0,0,1])
        
        sgj.plot3d()
        sgj.plot2d()
        
        Verbose.VERBOSITY=VERB_NONE

    def test_hulls(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        sg1=Spangler(nspangles=1000,sphash="Ring",n_equ=[1,0,5])
        sg1.populate_spangler(shape="ring",spangle_type=SPANGLE_GRANULAR,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)
        sg2=Spangler(nspangles=1000,sphash="Planet",n_equ=[0,0,1])
        sg2.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=1,seed=1,preset=True)
        sgj=Spangler(spanglers=[sg1,sg2])
        
        #Hulls of obsever
        cond,n_int,d_int=sgj.set_intersect(nvec=[1,0,0.1],center=[1,1,1]) #Each time a set intersect is executed the convex hulls are renewed
        sgj._calc_qhulls()
        
        fig,ax=plt.subplots()
        ax.scatter(sgj.data[cond].x_int,sgj.data[cond].y_int)
        f=convex_hull_plot_2d(sgj.qhulls["Planet"][0]["qhull"],ax)
        ax.axis("equal")

        f=convex_hull_plot_2d(sgj.qhulls["Ring"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][1]["qhull"],ax)

        #Hulls of light
        sgj.set_intersect([0,0,1]) #Each time a set intersect is executed the convex hulls are renewed
        sgj._calc_qhulls()
        fig,ax=plt.subplots()
        cond=sgj.data.visible
        ax.scatter(sgj.data[cond].x_int,sgj.data[cond].y_int)
        f=convex_hull_plot_2d(sgj.qhulls["Planet"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][1]["qhull"],ax)
        ax.axis("equal")
        
        Verbose.VERBOSITY=VERB_NONE

    def test_upint(self):
        plt.close("all")
        global sg
        
        Verbose.VERBOSITY=VERB_NONE
        
        """Shadow-test
        """
        nspangles=1000
        sps=[]
        sg=Spangler(nspangles=nspangles,sphash="Star",n_equ=[0,0,1],center_equ=[-7,0,0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_STELLAR,scale=3,seed=1,preset=1)
        sps+=[sg]
        nspangles=1000
        sg=Spangler(nspangles=nspangles,sphash="Planet",n_equ=[0,0,1])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=1,seed=1,preset=True)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Ring",n_equ=[1,0,2])
        sg.populate_spangler(shape="ring",spangle_type=SPANGLE_GRANULAR,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Moon",n_equ=[0,0,1],center_equ=[+4.0,0.0,0.0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_ATMOSPHERIC,scale=0.3,seed=1,preset=True)
        sps+=[sg]

        sg=Spangler(spanglers=sps)

        #"""
        sg.set_observer(nvec=sci.direction(30,0))
        sg.update_visibility_state()
        #""";

        #"""
        sg.set_luz(nvec=sci.direction(0,0))
        sg.update_illumination_state()
        #""";

        SHADOW_COLOR_LUZ=[90,0.2,1.0]
        sg.plot3d(center_at="Ring")
        sg.plot2d(center_at="Ring")
        sg._interact_plot2d(center_at="Ring",lon_luz=0)
        
        Verbose.VERBOSITY=VERB_NONE

    def test_anim(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE
        
        nspangles=500
        sps=[]
        sg=Spangler(nspangles=nspangles,sphash="Parent",n_equ=[0,0,1],center_equ=[-7,0,0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_STELLAR,scale=3,seed=1,preset=1)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Ring",n_equ=sci.direction(0,80))
        sg.populate_spangler(shape="ring",spangle_type=SPANGLE_GRANULAR,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Planet",n_equ=[0,0,1])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=1,seed=1,preset=True)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Moon",n_equ=[0,0,1],center_equ=[+3.0,0.0,0.0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_ATMOSPHERIC,scale=0.3,seed=1,preset=True)
        sps+=[sg]

        sg=Spangler(spanglers=sps)

        nobs=Plot.calc_flyby(normal=[0,0,1],start=0,stop=360,num=1,lat=30)
        nluz=Plot.calc_flyby(normal=[0,0,1],start=0,stop=360,num=5,lat=20)

        #anim,dirs=sg.animate_plot2d(nobs=nobs,nluz=nluz,interval=100,center_at="Ring",axis=False,not_plot=["Parent"])
        anim,dirs=sg._animate_plot2d(filename="/tmp/flyby-plot2d-luz.gif",nobs=nobs,nluz=nluz,interval=500,center_at="Ring",axis=False)
                
        Verbose.VERBOSITY=VERB_NONE

    def test_muluz(self):
        
        Verbose.VERBOSITY=VERB_NONE

        nspangles=100
        sps=[]

        sg=Spangler(nspangles=nspangles,sphash="Planet1",n_equ=[0,0,1],center_ecl=[0,0,0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=1,seed=1,preset=True)
        sps+=[sg]

        sg=Spangler(nspangles=nspangles,sphash="Moon1",n_equ=[0,0,1],center_ecl=[2,0,0])
        sg.populate_spangler(shape="sphere",spangle_type=SPANGLE_SOLID_ROCK,scale=0.5,seed=1,preset=True)
        sps+=[sg]

        sg=Spangler(spanglers=sps)

        sg.set_observer([1,1,1])
        sg.update_visibility_state()

        sphash="Planet1"
        sg.set_luz(nvec=[1,0,0],sphash=sphash)
        sg.update_illumination_state()

        #"""
        sphash="Moon1"
        sg.set_luz(nvec=[-2,1,0],sphash=sphash)
        sg.update_illumination_state()
        #"""

        sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
