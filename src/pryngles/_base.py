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
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module

# ## Module base
# 
# Goals of the module:
# - Future `base`
# - Basic Classes and functions of the Class.
# - It is based on existing infrastructure.

##HEADER
from pryngles import *
import unittest
#import rebound as rb
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

class Test(unittest.TestCase):
    from pryngles import _base

class Object(PrynglesCommon):
    """
    Global class of object
    """
    def _updateChilds(self,child=None):
        if 'childs' not in self.__dict__:
            self.childs=[]
        if child is not None:
            self.childs+=[child]
    def _updateParent(self,parent=None):
        if 'parent' not in self.__dict__:
            self.parent=parent
        elif parent is not None:
            self.parent=parent

def test_Object(self):
    obj=Object()
    obj._updateParent("parent")
    obj._updateChilds("child1")
    obj._updateChilds("child2")
    self.assertEqual([obj.parent],["parent"],True)
    print(obj.parent,obj.childs)
Test.test_Object=test_Object
unittest.main(argv=['first-arg-is-ignored'],exit=False)

def addStar(self,
                 hash=None,
                 center=None,
                 m=1,R=1,
                 Prot=0.1,
                 N=1000,
             star=None
            ):
    """
    TO DEVELOPERS:
        This method is for the System class.  It is placed here because when a new parameter
        is added to the initializer of the addObject (being Object Star, Planet, etc.) class
        then the parameters of this ruotine should also change.
    """
    if star is not None:
        self._addObject("star","Star",star)
    else:
        self.stars+=[Star(hash,center,m,R,Prot,N)]
    self._updateSystem()
    
class Star(Object):
    """
    Creates a star.
    """
    def __init__(self,
                 hash=None,
                 center=None,
                 m=1,R=1,
                 Prot=0.1,
                 N=1000,
                ):
        #List of arguments: hash,center,m,r,Prot,N
        #Hash of the object
        self.hash=hash
        self.type="Star"

        #Center of the system (object)
        self.center=center
        if self.center is not None:
            self.center._updateChilds(self)
            self._updateParent(self.center)
        
        #Basic common properties
        self.m=m #mass [cu]
        self.R=R #radius [cu]

        #Dynamic parameters
        self.Prot=Prot

        #Sampling parameters
        self.N=N #number of spangles
        
        #Update properties
        self.updateObject(**self.__dict__)

    def updateObject(self,**props):
        self.__dict__.update(props)
        self._updateChilds()
        self._updateParent()

def addPlanet(self,
                 hash=None,
                 center=None,
                 m=1e-3,R=0.1,
                 Prot=0.01,
                 N=1000,
             planet=None
            ):
    if planet is not None:
        self._addObject("planet","Planet",planet)
    else:
        self.planets+=[Planet(hash,center,m,R,Prot,N)]
    self._updateSystem()
    
class Planet(Object):
    """
    Creates a star.
    """
    def __init__(self,
                 hash=None,
                 center=None,
                 m=1e-3,R=0.1,
                 Prot=0.01,
                 N=1000,
                ):
        #List of arguments: hash,center,m,r,Prot,N
        #Hash of the object
        self.hash=hash
        self.type="Planet"

        #Center of the system (object)
        if center is None:
            raise ValueError(f"You must provide a valid object for the center.  {center} provided")
        self.center=center
        self.center._updateChilds(self)
        self._updateParent(self.center)
        
        #Basic common properties
        self.m=m #mass [cu]
        self.R=R #radius [cu]

        #Dynamic parameters
        self.Prot=Prot

        #Sampling parameters
        self.N=N #number of spangles
        
        #Update properties
        self.updateObject(**self.__dict__)
        
    def updateObject(self,**props):
        self.__dict__.update(props)
        self._updateChilds()
        self._updateParent()

def addRing(self,
                 hash=None,
                 center=None,
                 fi=1.5,fe=2.5,
                 N=1000,
             ring=None
            ):
    if ring is not None:
        self._addObject("ring","Ring",ring)
    else:
        self.rings+=[Ring(hash,center,fi,fe,N)]
    self._updateSystem()
    
class Ring(Object):
    """
    Creates a star.
    """
    def __init__(self,
                 hash=None,
                 center=None,
                 fi=1.5,fe=2.5,
                 N=1000,
                ):
        #List of arguments: hash,center,fi,fe,N
        #Hash of the object
        self.hash=hash
        self.type="Ring"

        #Center of the system (object)
        if center is None:
            raise ValueError(f"You must provide a valid object for the center.  {center} provided")
        self.center=center
        self.center._updateChilds(self)
        self._updateParent(self.center)

        #Basic common properties
        self.fi=fi
        self.fe=fe

        #Sampling parameters
        self.N=N #number of spangles
        
        #Update properties
        self.updateObject(**self.__dict__)
        
    def updateObject(self,**props):
        self.__dict__.update(props)
        self._updateChilds()
        self._updateParent()
        
        #Derivative properties
        self.ri=self.fi*self.center.R
        self.re=self.fe*self.center.R

System_doc="""
Creates a planetary system.

Examples:

    sys=System()
    sys=System(units=['msun','km','s'])
    sys=System(stars=dict(M=1,R=1))

Initialization (primary) attributes:

    units = ['msun','au','yr']: 
        Units used in calculations following the conventions and signs of rebound. List [3].
        
    stars = None: 
        Stars in the system.
    planets = None: 
        Planets in the system.
    rings = None: 
        Rings in the system.
           
        NOTE: Objects are provided either by description or as an object. Those attributes can be lists or 
              a single value.
    
    rebound = True:
        Set True if you want to simulte the orbit of objects (for instance if you want to calculate
        TTVs).  If False it will calculate orbits using Keplerian dynamics.

Other Public (secondary) attributes:

    N: Number of objects in the system.

Important private attributes:

    _sim: Rebound simulation object.
""";

class System(PrynglesCommon):
    
    def __init__(self,
                 units=['Msun','au','yr'],
                 stars=None,planets=None,rings=None,
                 rebound=True
                ):
        #Initialize rebound
        self.units=units

        #Add components
        self._list2Objects("stars","Star",stars)
        self._list2Objects("planets","Planet",planets)
        self._list2Objects("rings","Ring",rings)

        #Behavior
        self.rebound=True
        
        #Update system
        self._updateSystem()
    
    """
    #Freezer
    def addStar(self,star=None):
        self._addObject("star","Star",planet)
        self._updateSystem()
    
    def addPlanet(self,planet=None):
        self._addObject("planet","Planet",planet)
        self._updateSystem()
    
    def addRing(self,ring=None):
        self._addObject("ring","Ring",ring)
        self._updateSystem()
    """

    def _addObject(self,kind,objclass,obj):
        error=""
        try:
            obj.__dict__
            exec(f"self.{kind}s=self.{kind}s+[obj] if len(self.{kind}s)!=0 is not None else [obj]")
            cond=eval(f"obj.type!='{objclass}'")
            if cond:
                error=f"You cannot add a {objclass} with a {obj.type} object"
        except:
            #exec(f"self.{kind}s+=self.{kind}s+[{objclass}(**obj)] if self.{kind}s[0] is not None else [{objclass}(**{kind})]")
            exec(f"self.{kind}s+=self.{kind}s+[obj] if len(self.{kind}s)!=0 else [obj]")
        if error!="":
            exec(f"self.{kind}s.pop()")
            raise AssertionError(error)
        self._updateSystem()
        
    
    def _updateSystem(self):
        #Count components
        self.Nstars=self.Nplanets=self.Nrings=0
        for obj in "stars","planets","rings":
            exec(f"self.N{obj}=len(self.{obj})")
        self.N=self.Nstars+self.Nplanets+self.Nrings

    def _list2Objects(self,attr,kind,comps):
        if comps is None:
            exec(f"self.{attr}=[]")
        else:
            #Check if stars is a list
            try:comps[0]
            except:comps=[comps]
            exec(f"self.{attr}=[]")
            for comp in comps:
                try:
                    comp.__dict__
                    exec(f"self.{attr}+=[comp]")
                except:
                    exec(f"self.{attr}+=[{kind}(**comp)]")

System.addStar=addStar
System.addPlanet=addPlanet
System.addRing=addRing
System.__doc__=System_doc

