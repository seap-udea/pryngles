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
        
        Verbose.VERBOSITY=VERB_SIMPLE
        
        print("Basic definition:")
        sg=Spangler(nspangles=1,center_equ=[0,0,0],n_equ=[1,0,0])
        Misc.print_df(sg.data.head(1))
        return 

        print("\nCenter equ:")
        sg=Spangler(nspangles=3,center_equ=[0,0,1])
        Misc.print_df(sg.data.head(1))

        print("\nCenter ecl:")
        sg=Spangler(nspangles=3,center_ecl=[0,0,1])
        Misc.print_df(sg.data.head(1))

        print("\nRotation:")
        sg=Spangler(nspangles=3,w=30*Consts.deg,q0=40*Consts.deg)
        sg.set_positions(t=1)
        Misc.print_df(sg.data.head(1))

        print("\nJoin:")
        sg1=Spangler(sphash="Body 1",nspangles=3,w=40*Consts.deg)
        sg2=Spangler(sphash="Body 2",nspangles=3,w=30*Consts.deg)
        sg=Spangler(spanglers=[sg1,sg2])
        sg.set_positions(t=1)
        Misc.print_df(sg.data)

        Verbose.VERBOSITY=VERB_NONE

    def test_reset(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(nspangles=100)
        sg.reset_state()
        print_df(sg.data.head())

        Verbose.VERBOSITY=VERB_NONE

    def test_scale(self):

        Verbose.VERBOSITY=VERB_SIMPLE

        sg=Spangler(center_ecl=[1,1,1],center_equ=[1,1,1])
        print_df(sg.data)

        sg.set_scale(5)
        print_df(sg.data)

        Verbose.VERBOSITY=VERB_NONE

    def test_pop(self):
        Verbose.VERBOSITY=VERB_SIMPLE

        #No preset
        sg=Spangler(nspangles=850,n_equ=[1,0,0])
        sg.populate_spangler(geometry="ring",
                             spangle_type=GASEOUS_SPANGLE,
                             scale=2,seed=1,ri=0.2)
        sg.sample.plot()
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        print_df(sg.data.head(3))

        #Using preset
        sg=Spangler(nspangles=850)
        sg.populate_spangler(geometry="ring",
                             preset=True,
                             spangle_type=SOLID_SPANGLE,ri=0.2)
        sg.sample.plot()
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        print_df(sg.data.head(3))
    
        #Sphere
        sg=Spangler(nspangles=100)
        sg.populate_spangler(geometry="sphere",scale=3,seed=1)
        sg.sample.plot(spangled=dict(color='r',alpha=0.1))
        sg.sample.ax.set_title(f"N={sg.nspangles}")
        sg.sample.fig.tight_layout()
        print_df(sg.data.head(3))
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_plot3d(self):
        Verbose.VERBOSITY=VERB_SIMPLE

        #No preset
        sg=Spangler(nspangles=850,n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",preset=True,
                             spangle_type=GASEOUS_SPANGLE,
                             scale=2,ri=0.2)
        sg.data.illuminated=True
        cond=sg.data.x_ecl>0
        sg.data.loc[cond,"visible"]=True
        sg.data.loc[cond,"transmit"]=True
        sg.plot3d()
    
        #Sphere
        sg=Spangler(nspangles=100)
        sg.populate_spangler(geometry="sphere",preset=True,scale=3)
        sg.data.illuminated=True
        cond=sg.data.x_ecl>0
        sg.data.loc[cond,"visible"]=True
        sg.data.loc[cond,"transmit"]=True
        sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_setint(self):
        Verbose.VERBOSITY=VERB_SIMPLE

        #No preset
        sg=Spangler(nspangles=850,sphash="Ring")
        sg.populate_spangler(geometry="ring",preset=True,
                             spangle_type=GASEOUS_SPANGLE,
                             scale=2,ri=0.2)
        center=[1,1,-1]
        cond,n_int,d_int=sg.set_intersect(nvec=[1,0,1],center=center,sphash="Ring")
        
        sg.data.loc[cond,SPANGLER_COL_OBS]=sg.data.loc[cond,SPANGLER_COL_INT].values
        
        print_df(sg.data.head(1))
        print_df(sg.data.tail(1))

        from scipy.spatial import convex_hull_plot_2d
        fig,ax=plt.subplots()

        f=convex_hull_plot_2d(sg.qhulls["Ring"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sg.qhulls["Ring"][1]["qhull"],ax)
        
        #Remove points corresponding to qhull
        for l in fig.axes[0].get_children():
            if type(l) is Line2D:
                plt.setp(l,ms=0,zorder=100)
    
        ax.scatter(sg.data.x_int,sg.data.y_int,color='c',s=65,fc="None",alpha=0.2,zorder=100)        
        ax.axis("equal")

        #Plot 3d
        sg.plot3d(coords="int")
        plane=sg.qhulls["Ring"][0]["plane"]
        plane.plot_plane(ax=sg.ax3d,color='r',alpha=0.5)
        
        #Hulls
        print(sg.qhulls)
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_setobs(self):
        Verbose.VERBOSITY=VERB_SIMPLE

        #No preset
        sg=Spangler(nspangles=850,n_equ=[1,1,2])
        sg.populate_spangler(geometry="ring",preset=True,
                             spangle_type=GRANULAR_SPANGLE,
                             scale=2,ri=0.2)
        sg.set_observer(nvec=[1,1,1],center=[0,0,1])
        sg.data.illuminated=True
        sg.data.loc[sg.data.cos_obs>0,"visible"]=True
        sg.plot3d()
    
        #Sphere
        sg=Spangler(nspangles=100)
        sg.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        sg.set_observer(nvec=[1,1,1])
        sg.data.illuminated=True
        sg.data.loc[sg.data.cos_obs>0,"visible"]=True
        sg.plot3d()
        
        #Test center
        sg=Spangler(nspangles=100,center_ecl=[0,0,3])
        sg.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        center=[0,0,-1]
        sg.set_intersect(nvec=[0,0,1],center=center)
        sg.data.illuminated=True
        sg.data.loc[sg.data.cos_int>0,"visible"]=True
        sg.plot3d(coords="int")
         
        Verbose.VERBOSITY=VERB_NONE
        
    def test_setluz(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        #Sphere
        sg=Spangler(nspangles=500)
        sg.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
 
        sg.set_observer([0,0,1])
        sg.data.loc[sg.data.cos_obs>0,"visible"]=True

        sg.set_luz([0,0,1])
        sg.data.loc[sg.data.cos_luz>0,"illuminated"]=True

        sg.plot3d()
        
        #Ring
        sg=Spangler(nspangles=850)
        sg.populate_spangler(geometry="ring",preset=True,
                             spangle_type=GRANULAR_SPANGLE,
                             scale=2,ri=0.2)

        sg.set_observer([1,1,+1])
        sg.data.visible=True

        sg.set_luz([0,0,+1])
        sg.data.illuminated=True
        sg.data.loc[(sg.data.cos_luz*sg.data.cos_obs)<0,"transmit"]=True
        
        sg.plot3d()
        
        #Two spheres
        sg1=Spangler(sphash="Star 1",nspangles=100,center_equ=[-5,0,0])
        sg1.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        
        sg2=Spangler(sphash="Star 2",nspangles=100,center_equ=[+5,0,0])
        sg2.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        
        sg=Spangler(spanglers=[sg1,sg2])

        sg.set_observer([0,1,0])
        sg.data.loc[sg.data.cos_obs>0,"visible"]=True

        sg.set_luz(nvec=[1,0,0],sphash="Star 1")
        sg.set_luz(nvec=[0,0,1],sphash="Star 2")
        sg.data.loc[sg.data.cos_luz>0,"illuminated"]=True

        sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_simplevis(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        plt.close("all")
        #Ring with semitransparent spangle: all illuminated, all visible, no transmission
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",ri=0.3,spangle_type=GRANULAR_SPANGLE,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([1,1,-1])
        sg.update_simple_state()
        sg.plot3d()

        #Ring with semitransparent spangle: all illuminated, all visible, no transmission
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(geometry="ring",ri=0.3,spangle_type=GRANULAR_SPANGLE,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([-1,-1,-1])
        sg.update_simple_state()
        sg.plot3d()
        
        #Sphere with solid spangle: only illuminated 
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,preset=True,scale=3)
        sg.set_observer([1,0,1])
        sg.set_luz([0,0,1])
        sg.update_simple_state()
        sg.plot3d()
        
        #Sphere with stellar spangle: all illuminated, not all visible
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",spangle_type=STELLAR_SPANGLE,preset=True,scale=3)
        sg.set_observer([1,0,1])
        sg.set_luz([0,0,1])
        sg.update_simple_state()
        sg.plot3d()

        #Sphere with semitransparent spangle: all illuminated, all visible
        sg=Spangler(nspangles=100,n_equ=[1,1,1])
        sg.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        sg.set_observer([0,0,1])
        sg.set_luz([1,0,0])
        sg.update_simple_state()
        sg.plot3d()

        #Two spheres
        sg1=Spangler(sphash="Star 1",nspangles=100,center_equ=[-5,0,0])
        sg1.populate_spangler(geometry="sphere",spangle_type=GASEOUS_SPANGLE,preset=True,scale=3)
        
        sg2=Spangler(sphash="Star 2",nspangles=100,center_equ=[+5,0,0])
        sg2.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,preset=True,scale=3)
        
        sg=Spangler(spanglers=[sg1,sg2])

        sg.set_observer([0,1,0])
        sg.set_luz(nvec=[1,0,0],sphash="Star 1")
        sg.set_luz(nvec=[0,0,1],sphash="Star 2")
        sg.update_simple_state()
        
        sg.plot3d()
        
        Verbose.VERBOSITY=VERB_NONE
        
    def test_plotobs(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE
        
        plt.close("all")
        sg=Spangler(nspangles=2500,sphash="123",n_equ=[1,1,1],center_ecl=[0,0,2])
        sg.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=2,seed=1,preset=True)
        sg.set_observer(nvec=[1,0,0])
        sg.set_luz(nvec=[1,1,1])
        sg.update_simple_state()
        sg.plot_obs()
        
        Verbose.VERBOSITY=VERB_NONE

    def test_join(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        sg1=Spangler(nspangles=1000,sphash="Ring",n_equ=[1,0,5])
        sg1.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)

        sg2=Spangler(nspangles=1000,sphash="Planet",n_equ=[0,0,1])
        sg2.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=1,seed=1,preset=True)

        sgj=Spangler(spanglers=[sg1,sg2])
        
        sgj.set_observer([1,0,0.1])
        sgj.set_luz([0,0,1])
        sgj.update_simple_state()
        
        sgj.plot3d()
        sgj.plot_obs()
        
        
        Verbose.VERBOSITY=VERB_NONE

    def test_hulls(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE

        sg1=Spangler(nspangles=1000,sphash="Ring",n_equ=[1,0,5])
        sg1.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)
        sg2=Spangler(nspangles=1000,sphash="Planet",n_equ=[0,0,1])
        sg2.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=1,seed=1,preset=True)
        sgj=Spangler(spanglers=[sg1,sg2])
        
        #Hulls of obsever
        cond,n_int,d_int=sgj.set_intersect(nvec=[1,0,0.1],center=[1,1,1]) #Each time a set intersect is executed the convex hulls are renewed
        
        fig,ax=plt.subplots()
        ax.scatter(sgj.data[cond].x_int,sgj.data[cond].y_int)
        f=convex_hull_plot_2d(sgj.qhulls["Planet"][0]["qhull"],ax)
        ax.axis("equal")

        f=convex_hull_plot_2d(sgj.qhulls["Ring"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][1]["qhull"],ax)

        #Hulls of light
        sgj.set_intersect([0,0,1]) #Each time a set intersect is executed the convex hulls are renewed
        fig,ax=plt.subplots()
        cond=sgj.data.visible
        ax.scatter(sgj.data[cond].x_int,sgj.data[cond].y_int)
        f=convex_hull_plot_2d(sgj.qhulls["Planet"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][0]["qhull"],ax)
        f=convex_hull_plot_2d(sgj.qhulls["Ring"][1]["qhull"],ax)
        ax.axis("equal")
        
        Verbose.VERBOSITY=VERB_NONE

    def test_plotint(self):
        
        Verbose.VERBOSITY=VERB_SIMPLE
        
        plt.close("all")

        nspangles=100
        sps=[]
        sg=Spangler(nspangles=nspangles,sphash="Parent",n_equ=[0,0,1],center_equ=[-5,0,0])
        sg.populate_spangler(geometry="sphere",spangle_type=STELLAR_SPANGLE,scale=2,seed=1,preset=True)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Planet",n_equ=[0,0,1])
        sg.populate_spangler(geometry="sphere",spangle_type=SOLID_SPANGLE,scale=1,seed=1,preset=True)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Ring",n_equ=[1,0,0])
        sg.populate_spangler(geometry="ring",spangle_type=GRANULAR_SPANGLE,scale=2.5,seed=1,ri=1.5/2.5,boundary=0)
        sps+=[sg]
        sg=Spangler(nspangles=nspangles,sphash="Moon",n_equ=[0,0,1],center_equ=[+3.0,0.0,0.0])
        sg.populate_spangler(geometry="sphere",spangle_type=LIQUID_SPANGLE,scale=0.3,seed=1,preset=True)
        sps+=[sg]
        sg=Spangler(spanglers=sps)
        n_obs=sci.cartesian([1,45*Consts.deg,10*Consts.deg])
        cond=sg.set_observer(nvec=n_obs,center=[0,0,0])
        
        #Plot intersect
        sg.reset_state()
        sg.update_simple_state()
        sg.update_visibility_state()
        fig=sg._plot_intersect(prop="visible")
        
        #View intersect
        sg._view_intersect(lon=45,lat=30)
        
        #Interact intersect
        sg._interact_intersect()
        
        #Animation
        sg._animate_intersect(filename="/tmp/test.mp4",lat=5,num=5)
        
        Verbose.VERBOSITY=VERB_NONE


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
