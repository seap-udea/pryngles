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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External required packages
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

from pryngles import *

import rebound as rb
from tqdm import tqdm


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class System
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class System(PrynglesCommon):
    """
    High-level interface for building spangled planetary systems.
    Integrates Rebound simulations with the Pryngles spangling pipeline to generate discrete surface facets,
    manage bodies, compute photometric states, and evolve dynamics.

    Parameters
    ----------
    filename : `str`
        Path to a saved system file to load. If provided, units and resetable
        parameters are ignored | **Default:** None
    units : `list` 
        Three-element list of str specifying units for length, mass, and time
        as recognized by Rebound (order: length, mass, time) | **Default:** ['au', 'msun', 'yr2pi']
    resetable : `bool`
        Whether the system state can be snapshotted and reset to initial | **Default:** False

    Attributes
    ----------
    sim : `rebound.Simulation`
        Underlying Rebound simulation instance (None until initialized).
    ul, um, ut : `float`
        Conversion factors from internal units to SI (meters, kilograms, seconds).
    G : `float`
        Gravitational constant in chosen units.
    bodies : `odict`
        Ordered dictionary of :any:`body.Body` instances in the system, keyed by name.
    nbodies : `int`
        Number of bodies in `bodies`.
    nparticles : `int`
        Number of particles in the Rebound simulation.
    n_obs : `np.array(3)`
        Normalized observer direction vector :math:`\hat{n}` in the system's coordinate frame.
    sg : :any:`spangler.Spangler`
        Spangler instance covering all system facets.
    spangle_scatterers : `dict`
        Mapping from :any:`spangler.Spangler` type attribute to :any:`scatterer.Scatterer` classes and options.
        This means that for spangles of the type :any:`consts.SPANGLE_ATMOSPHERIC` Pryngles will 
        instantiate an object of the class :any:`scatterer.LambertianGrayAtmosphere` to compute their scattering.

    Examples
    --------
    >>> # Create a new system with default units
    >>> sys = System(units=['au','msun','yr2pi'])
    >>> 
    >>> # Now you can add bodies to the system
    >>> star = sys.add('Star', m=1.0)
    >>> planet = sys.add('Planet', parent=star, m=1e-3, a=1.0, e=0.1)
    >>>
    >>> # Initialize the simulation
    >>> sys.initialize_simulation()
    >>>
    >>> # Spangle the bodies to compute photometry
    >>> sys.spangle_system()
    >>> 
    >>> # You can visualize the system
    >>> sys.sg.plot2d(include = ['Star', 'Planet'])

    .. image:: images/system_init.png
        :align: center
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Bassic methods
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def __init__(self,
                 filename=None,
                 units=['au','msun','yr2pi'],
                 resetable=False
                ):
        
        if filename:
            self.load_from(filename)
            return
        
        #Rebound simulation
        self.sim=None
        self._simulated=False
        
        #Attributes by default
        
        #List of bodies in the system
        self.bodies=odict()
        
        #Root of the tree of bodies
        self.root=None
        
        #Center of the light-source in the system
        self.source=None
        self.center_root=np.array([0,0,0])
        
        #Orbital configuration
        self.orbital_configuration=None
        
        #Observer properties
        self.n_obs=np.array([0,0,1])
        self.alpha_obs=0  
        self.center_obs=None
        
        #Check if spangled
        self._spangled=False
        
        #Check if observer has been set
        self._observer_set=False
        self._luz_set=False
        
        #Initialize spangler object
        self.sg=None
        
        #Is the system resetable?
        self._resetable=resetable
        if self._resetable:
            #Create temporary file
            self._snap_file_name = "/tmp/pryngles-system.pkl"
        
        #Update rebound units
        self.update_units(units)
        
        #By default spangle scatterers
        """
        This is the list of the class of scatterers used to calculate the scattering in
        different types of spangles.
        
        The structure of this dictionary is:
        
            key: integer (or enumerator):
                This is the column spangle_type in the spangler DataFrame.
                
            item: tuple (2):
                Component 1: 
                    Class of scatterer.
                Component 2: 
                    Dictionary mapping the initialization properties of the scatterer to columns
                    in the spangler DataFrame.
        
        Example of item:
        
            SPANGLE_ATMOSPHERIC:(LambertianGrayAtmosphere,dict(AS="albedo_gray_spherical"))
            
                This means that for spangles of the type SPANGLE_ATMOSPHERIC Pryngles will 
                instantiate an object of the class LambertianGrayAtmosphere.  This class have a
                single parameter, the spherical albedo AS.  The dictionary means that when 
                instantiating the object the column "albedo_gray_spherical" will be used to 
                initialize the object.                
        """
        self.spangle_scatterers={
            SPANGLE_ATMOSPHERIC:(LambertianGrayAtmosphere,dict(AS="albedo_gray_spherical")),
            SPANGLE_GRANULAR:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_LIQUID:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_SOLID_ICE:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_SOLID_ROCK:(LambertianGraySurface,dict(AL="albedo_gray_normal")),
            SPANGLE_GASEOUS:(BlackBodySurface,dict()),
            SPANGLE_STELLAR:(BlackBodySurface,dict()),
        }
        
    def update_units(self,units):
        """
        Set or change the internal unit conversion factors.

        Parameters
        ----------
        units :  `list`
            Three-element string of Rebound-recognized units: [length, mass, time]. 
            You can see the defined ones in :any:`consts.REBOUND_ORBITAL_PROPERTIES`

        Raises
        ------
        ValueError
            If any unit string is not recognized by Rebound.
        """
        #Check units
        if units[0] not in rb.units.lengths_SI:
            raise ValueError(f"Length unit provided '{units[0]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.lengths_SI.keys())}")
        if units[1] not in rb.units.masses_SI:
            raise ValueError(f"Mass unit provided '{units[1]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.masses_SI.keys())}")
        if units[2] not in rb.units.times_SI:
            raise ValueError(f"Time unit provided '{units[2]}' is not recognized by Rebound.  Use one of these: {tuple(rb.units.times_SI.keys())}")
        
        #Units        
        self.units=units
        self._ul,self._um,self._ut=self.units
        #self.sim.units=self.units
        
        #Canonical units of the system
        self.ul=rb.units.convert_length(1,self._ul,"m")
        self.um=rb.units.convert_mass(1,self._um,"kg")

        #Compute the units of time
        sim=rb.Simulation()
        sim.units=self.units
        self.G=sim.G
        self.ut=np.sqrt(self.G*self.ul**3/(self.um*GSI))
        
        #Update system
        self._update_system()
        
    def _get_source(self,body):
        """Get the source of light (stellar body) in the center of a body
        """
        if (body.parent is None) or (body.kind == "Star"):
            return body

        elif body.parent.kind == "Star":
            return body.parent

        else:
            return self._get_source(body.parent)

    def _update_system(self):
        """Update system properties
        """
        self.nbodies=len(self.bodies)
        if self._simulated:
            self.nparticles=len(self.sim.particles)
        
    def _is_spangled(self):
        """Check if system is spangled
        """
        return True if self.sg else False
    
    def reset_state(self):
        """
        Restore Spangler visibility and illumination to initial.
        Calls :meth:`Spangler.reset_state` and clears observer/light flags.
        """
        self.sg.reset_state()
        self._observer_set=False
        self._luz_set=False

    def save_to(self,filename):
        """Save system from file
        
        Parameters
        ----------
        filename : str
            Base filename (without extension) for saving.

        Returns
        --------------
        :
            filename : str
             Rebound state is saved as filename + '.rbin'.
        """
        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Save rebound state
            verbose(VERB_SIMPLE,"Saving rebound simulation")
            self.sim.save(rb_filename)

            #Since rebound have ctypes it cannot be pickled
            del self.sim
            self._simulated=True

        #Pickle system
        PrynglesCommon.save_to(self,filename)

        if self._simulated:
            #Load again rebound
            self.sim=rb.Simulation(rb_filename)

    def load_from(self,filename):
        """Load system from filename
                
        Parameters
        ----------
        filename : str
            Base filename (without extension) used for pickle and Rebound state.

        Returns
        ------------
        :
            sim : ``rebound.Simulation``
             :data:`~ system.System.sim` attribute loaded from file
        """
        #Load system
        self=PrynglesCommon.load_from(self,filename)

        if self._simulated:
            #Rebound file
            rb_filename=filename+".rbin"

            #Load rebound
            verbose(VERB_SIMPLE,"Loading rebound simulation")
            self.sim=rb.Simulation(rb_filename)
        
    def status(self):
        """
        Print a summary of the current system and particles.

        Note
        ---------
        If not yet initialized, prompts for initialization.
        """
        if self._simulated:
            print(f"System with {self.nbodies} bodies and {self.nparticles} particles (rings and disk are not particles)")
            self.sim.status()
        else:
            print(f"Simulation for this system has not been yet initialized. Use System.initialize_simulation()")


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file scatterer
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def update_scatterers(self):
        """
        Instantiate scatterer objects for all spangles lacking one.

        Raises
        ------
        AssertionError
            If spangles have not yet been generated via :data:`~ system.System.spangle_system()`.
        """
        if not self._spangled:
            raise AssertionError("You need to spangle the system before updating the scatterers.")
        
        #Update scatterer only for the non-assigned one
        cond=(self.data.scatterer=="")
        
        for index in self.data[cond].index:
            #Get spangle
            spangle=self.data.loc[index]
            #Get spangle sype
            spangle_type=spangle["spangle_type"]
            #Get scatterer class and options description
            spangle_scatterer,spangle_options=self.spangle_scatterers[spangle_type]
            #Build options of scatterers from options description
            scatterer_options={**dict(zip(spangle_options.keys(),spangle[list(spangle_options.values())]))}
            #Instantiate object of scatterer and save hash into DataFrame
            self.data.loc[index,"scatterer"] = spangle_scatterer(**scatterer_options) #.hash To not save hash but the Scatterer object
    
    def update_albedos(self):
        """ 
        Compute directional-dependent Lambertian albedo per spangle. 
        It implements :data:`~ scatterer.Scatterer.get_albedo()` method. See our :doc:`scatterer` for the theory behind our models

        Note
        ------
        Applies only to non-stellar spangles, i.e., with no :data:`~ consts.SPANGLE_STELLAR` attribute
        """

        # Creating lambertian_albedo column for directional-dependent albedo for surface/atmosphere scattering
        self.data['lambertian_albedo'] = 0.0

        # Only planetary surfaces are taken into account
        cond = (self.data['spangle_type'] != 6) & (self.data['cos_luz'] >= 0)
    
        # Initialize Lambertian Albedo columns into de DataFrame
        self.data.loc[cond, "lambertian_albedo"] = self.data.loc[cond].apply(
            lambda sp: sp["scatterer"].get_albedo(
                eta = sp["cos_obs"], zeta = 0, delta = 0, lamb = 0.55), axis = 1)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Tested methods from module file system
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def add(self,kind="Star",parent=None,**props):
        """
        Add a new Body to the System and auto-spangle it.

        Parameters
        ----------------
        kind : `str`
            Type of body to add from :any:`consts.BODY_KINDS` | **Defaults** = `Star`
        parent : Body or None
            Parent body (orbital central object); None for the root star.
        **props : dict
            Properties forwarded to the Body constructor (mass, radius, orbital elements). It came from :any:`consts.BODY_DEFAULTS`

        Returns
        -------
        :
            output : :any:`body.Body`
             The newly created Body instance.

        Raises
        ------
        AssertionError
            If kind is None.
        ValueError
            If kind not in :any:`consts.BODY_KINDS` or name already exists.    
        """
        if kind is None:
            raise AssertionError("You must provide a valid object kind (Star, Planet, Ring).")
    
        if kind not in BODY_KINDS:
            raise ValueError(f"Object kind '{kind}' is not recognized.")
    
        #Legacy
        if 'primary' in props:
            parent=props["primary"]
        if kind=="Observer":
            parent=self.root
        
        #Default parameters
        if self.root:
            if (kind!="Star") and (parent is None):
                parent=self.root
                
            if kind=="Planet":
                if "m" not in props:
                    props["m"]=0.1*parent.m
                if "radius" not in props:
                    props["radius"]=0.5*parent.radius
                if "a" not in props:
                    if "a" in parent.__dict__:
                        props["a"]=0.5*parent.a
                    
        #Create body
        props.update(dict(name_by_kind=True))
        self.__body=eval(f"{kind}(parent=parent,**props)")
        
        if self.__body.name in self.bodies:
            raise ValueError(f"An object with name '{self.__body.name}' has been already added.")
        
        #If we have a root object and no parent has been provided
        if not parent:
            if self.root:
                raise ValueError(f"A root object alread exist in the system ({self.root.name}) and you do not provided a parent body for {self.__body.name}.")
            else:
                self.root=self.__body
                verbose(VERB_SIMPLE,f"Setting the root object as {self.root.name}")
            
        self.bodies[self.__body.name]=self.__body
        
        #Create the shined body tree
        self.__body.shined=[]
    
        #Update system
        self._update_system()
        
        #Set the source of the object
        if self.__body.source:
            #Check that the source is a body
            if not isinstance(self.__body.source,Body):
                raise ValueError(f"The source of body must be an actual Body.")        
            #Check that the source is among the bodies
            if self.__body.source.name not in self.bodies:
                raise ValueError(f"The source of body {self.__body.name} is not among system bodies {list(self.bodies.keys())}.")
            #Check that the source be a star
            if self.__body.source.kind!="Star":
                raise ValueError(f"The source of body {self.__body.name} must be a Star.  You set {self.__body.source.name} which is a {self.__body.source.kind}.")
            self.__body.source.shined+=[self.__body.name]
        else:
            if self.__body.kind=="Star":
                self.__body.source=self.__body
            elif self.__body.parent==self.root:
                self.__body.source=self.root
            else:
                self.__body.source=self.__body.parent.source
        self.__body.source.shined+=[self.__body.name]
        
        verbose(VERB_SIMPLE,f"Object '{kind}' with name '{self.__body.name}' has been added.")

        # if kind == 'Ring':
        #     self.__body.m = -1

        return self.__body
    
    
    def initialize_simulation(self,orbital_tree=None,**rebound_options):
        """
        Build and initialize the Rebound simulation from the orbital tree scheme.

        Note
        ----------
        See the :doc:`orbit` for the convenction to construc an orbital tree input

        Parameters
        ----------
        orbital_tree : `list`
            Hierarchical tree of bodies for N-body initialization; see :any:`orbit.OrbitUtil.build_tree`
        
        Returns
        -------
        :
            output : :any:`orbit.Orbit`
             The Orbit object representing the hierarchical system.

        Raises
        ------
        AssertionError
            If a Body is in System but not in orbital tree.
        """

        # See Orbit Module for build_tree() method
        if orbital_tree is None: self.orbital_tree = OrbitUtil.build_tree(self.root)

        else: self.orbital_tree = orbital_tree
        
        #Set the rebound hash of all bodies
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                body.rbhash=body.parent.name
            else:
                body.rbhash=body.name
        
        #Check that all bodies in system is in the orbital tree
        bodies=list(Misc.flatten(self.orbital_tree))
        for name,body in self.bodies.items():
            if body.kind=="Ring":
                continue
            if body not in bodies:
                raise AssertionError(f"Body '{name}' is in System but not in orbital tree.")
            
        #Build hierarchical N-body system
        orbit,pelements=OrbitUtil.build_system(self.orbital_tree,self.units)
        orbit.calculate_orbit()
        
        #Initialize simulation
        self.sim=rb.Simulation()
        self.sim.units=self.units
        
        #Add particles to simulation
        orbit.sim.move_to_com()
        for i,p in enumerate(orbit.sim.particles):
            self.sim.add(
                hash=bodies[i].name,
                m=bodies[i].m,
                x=p.x,y=p.y,z=p.z,
                vx=p.vx,vy=p.vy,vz=p.vz
            )
        self.sim.orbit=orbit
        self._simulated=True
        self._update_system()
        
        return orbit
    
    def remove(self,name):
        """
        Remove a body (and children objects) from the system.

        Parameters
        ----------
        name : str
            Name (hash) of the body to remove.

        Raises
        ------
        ValueError
            If no body with the given name exists.
        """
        
        if name in self.bodies:
            verbose(VERB_SIMPLE,f"Removing object {name} from system")
    
            obj=self.bodies[name]
    
            #Get the list of child hashes before removing (it changes during for)
            child_hashes=list(obj.childs.keys())
            
            #Remove child objects
            for child_hash in child_hashes:
                if child_hash in self.bodies:
                    self.remove(child_hash)
                    
            #Remove object from Rebound simulation
            if obj.kind != "Ring":
                if self._simulated:
                    if self.nparticles:
                        verbose(VERB_SIMPLE,f"Removing particle {name} from simulation")
                        self.sim.remove(hash=name)
            
            #Remove object from childs of its parent
            if obj.parent:
                del obj.parent.childs[name]
            
            #Remove object from bodies
            del self.bodies[name]
    
            #Update system
            self._update_system()
        else:
            raise ValueError(f"No object with hash '{name}' in the system")
    
    def spangle_system(self):
        """
        Generate :any:`spangler.Spangler` instances for all bodies and join them.

        Attributes
        ------------
        sg : :any:`spangler.Spangler`
            It contains the :data:`~ spangler.Spangler` object in wich we sample and discretize the surface of all bodies defined in system in order to compute light-matter interactions
        data : `pd.DataFrame`
            Contains :any:`consts.SPANGLER_COLUMNS` and sets the  default observer/light states.

        Raises
        ------
        AssertionError
            If simulation not yet initialized.
        """
        if not self._simulated:
            raise AssertionError("Before spangling the system you must initialize the simulation: System.initialize_simulation().")
        
        self._spanglers=dict()
        
        #Add spangles
        for name,body in self.bodies.items():
            
            verbose(VERB_SIMPLE,f"Spangling body '{name}' (kind '{body.kind}')")
            body.spangle_body()
    
            #Center object around its position according to rebound
            body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
            body.sg.set_positions(center_ecl=body.center_ecl)
            self._spanglers[name]=body.sg
            
        #Set the center of the source of light for each body
        for name,body in self.bodies.items():
            body.center_source=body.source.center_ecl
            if body==self.root:
                self.center_root=body.source.center_ecl
                
        #Join spanglers
        self.sg=Spangler(spanglers=list(self._spanglers.values()))
    
        #An usefule alias
        self.data=self.sg.data
        
        #Set default observer
        self.update_perspective(n_obs=self.n_obs,alpha_obs=self.alpha_obs)
        
        #Save state of the system
        if self._resetable:
            self.save_to(self._snap_file_name)
        
        #Already spangled
        self._spangled=True

        #Update Scatterers for Albedo Computing
        self.update_scatterers()
    
    def _set_observer(self,nvec=[0,0,1],alpha=0,center=None):
        """Set the position of the observer
        """
        #Only set observer if it is spangled
        if self._is_spangled():
            
            #At changing the observer, reset state
            self.sg.reset_state()
            
            #Set observer
            self.sg.set_observer(nvec=nvec,alpha=alpha,center=center)
            
            #Update areas of the spangles
            
            #Update system properties
            self.d_obs=self.sg.d_obs
            self.n_obs=self.sg.n_obs.copy()
            self.rqf_obs=self.sg.rqf_obs.copy()
            self.alpha_obs=self.sg.alpha_obs
            self.center_obs=self.sg.center_obs
        
            #Update visibility
            self.sg.update_visibility_state()
            
            #Check that observer has been set
            self._observer_set=True
            
        else:
            raise AssertionError("You must first spangle system before setting observer direction.")
            
    def _set_luz_recursive(self,name,nluz,verbosity=VERB_SIMPLE):
        """Set light source for body and 
        """
        body=self.bodies[name]
        verbose(verbosity,f"Illuminating body {name}, with nluz = {nluz} and center = {body.center_source}")
        self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
        if body.childs:
            verbose(verbosity,f"\tObject {name} has childs!")
            for child_name in body.childs:
                verbose(verbosity,f"\tCalling recursively set_luz for {child_name}")
                self._set_luz_recursive(child_name,nluz)
        else:
            verbose(verbosity,f"\tObject {name} has no childs!")
                
    def _set_luz(self,verbosity=VERB_SIMPLE):
        """Set illumination in the system.
        
        Update:
            States: illuminated, shadow, hidden_by_luz
        """
        if self._is_spangled():
            
            if not self._observer_set:
                raise AssertionError("You must first set observer before setting light.")
            
            self.bodies_illuminated=[]
            for name,body in self.bodies.items():
              
                if body.kind == "Star":
                    verbose(verbosity,f"Body {body.name} is a star... skipping")
                    continue
                    
                if body.parent.name in self.bodies_illuminated:
                    verbose(verbosity,f"Parent body of {name}, {body.parent.name}, has been already illuminated")
                    continue
                            
                #Get center of body
                center=body.center_ecl
                        
                #Get source and center
                verbose(verbosity,f"Calculating illumination for '{name}' coming from '{body.source.name}' @ {body.center_source}")            
                nluz=body.center_source-center
                        
                if body.kind == "Ring" and body.parent.kind == "Star":
                    verbose(verbosity,f"Parent body of ring, {body.parent.name} is a star. All spangles will be illuminated")
                    self.sg.set_luz(nvec=nluz,center=body.center_source,name=name)
                    cond=(self.sg.data.name==name)
                    self.sg.data.loc[cond,"unset"]=False
                    self.sg.data.loc[cond,"illuminated"]=True
                    self.sg.data.loc[cond,"shadow"]=False
                else:                
                    verbose(verbosity,f"Illuminating body {name} and all its childs")
                    self._set_luz_recursive(name,nluz,verbosity)
                    self.sg.update_illumination_state(included=body.source.shined)
                    self.bodies_illuminated+=[name]
                    
            self._luz_set=True
        else:
            raise AssertionError("You must first spangle system before setting light.")
            
    def update_perspective(self,n_obs=None,alpha_obs=0,center_obs=None):
        """
        Set new observer and illumination perspectives for direction values and recompute photometry.

        Parameters
        ----------
        n_obs : `array_like`
            Observer direction vector.
        alpha_obs : `float`
            Observer roll angle (radians).
        center_obs : `array_like`
            Observer position vector.
        
        Warning
        ----------
        This function may have slow performance
        """
        if n_obs is not None:
            
            #Update observing conditions
            self.n_obs,one=spy.unorm(n_obs)
            self.alpha_obs=alpha_obs
            self.center_obs=center_obs
    
        #Set observer
        self._set_observer(nvec=self.n_obs,alpha=self.alpha_obs,center=center_obs)
        self._set_luz()

        #Set Photomety Values
        # self.update_StellarFlux()
        # self.update_DiffuseReflection()
        
    
    def update_body(self,body,**props):
        """
        Update a body's properties before spangling.

        Parameters
        ----------
        body : `str` or :data:`~ body.Body`
            Identifier or instance of the body to update.
        **props : dict
            Properties to update (mass, radius, orbital elements).

        Raises
        ------
        AssertionError
            If system already spangled or body not found.
        """
        #Update spangling?
        if self._is_spangled():
            raise AssertionError("After spangling you cannot update the properties of the bodies.  Please rebuild the system")
    
        #Update body properties
        if isinstance(body,Body):
            body.update_body(**props)
        elif body in self.bodies:
            body=self.bodies[body]
            lkind=body.kind.lower()
            exec(f"body.update_{lkind}()")
        else:
            raise AssertionError("You are trying to update a body ({body}) which is not in the system")
            
        #Check if among props there is any property related to position
        if any(k in props for k in REBOUND_ORBITAL_PROPERTIES):
            raise ValueError(f"You cannot update an orbital property {props} without compromising the full simulation. Rebuild the system from scratch.")
    
    def reset(self):
        """
        Reload system state from the last snapshot.
        Works only if ``resetable = True`` on creation.
        """
        if self._resetable:
            self.load_from(self._snap_file_name)
            pass
        else:
            print("System is not resetable. Use resetable = True when defining the System or when you spangle it.")
    
    def integrate(self,*args,**kwargs):
        """
        Advance N-body simulation by a given time step, update center of each body and set positions of the Spangles.

        Parameters
        ----------
        t : `float`
            Time for integrate to it

        Raises
        ------
        AssertionError
            If system not yet spangled.
        
        Note
        ---------
        It doesn´t update the Visibility and Illumination states of the Spangles.
        """
        #Time of integration
        t=args[0]
        verbose(VERB_SIMPLE,"Integrating up to {t}")
        
        if self._spangled:
            
            #Integrate
            self.sim.integrate(*args,**kwargs)
            self.sim.move_to_com()
        
            #Update positions
            for name,body in self.bodies.items():
                
                #Position of the body according
                body.center_ecl=np.array(self.sim.particles[body.rbhash].xyz)
    
                verbose(VERB_VERIFY,f"Updating center of body {name} @ {body.center_ecl}")
                cond=self.sg.data.name==name
                self.sg.data.loc[cond,"center_ecl"]=pd.Series([list(body.center_ecl)]*sum(cond),dtype=object).values
    
            #Update positions
            self.sg.set_positions()
            
        else:
            raise AssertionError("You must first spangle system before setting positions.")

    def integrate_perspective(self, t, n_obs = None, alpha_obs = 0, center_obs = None):
        """
        Rebound Integration of N-Body System and update center of each body and set positions of the Spangles. 
        Also update Visibility and Illumination states of the Spangles


        Parameters
        ----------
        t : float
            Target time for integration.
        n_obs : array_like or None
            Observer direction.
        alpha_obs : float
            Observer roll angle.
        center_obs : array_like or None
            Observer position.

        Warning
        --------------
        This function may have slow performance
        """
        #Integration
        self.integrate(t)

        #Updating Spangles Perspective
        self.update_perspective(n_obs = None, alpha_obs = 0, center_obs = None)
    
    def ensamble_system(self,lamb=0,beta=0,**physics):
        """
        Legacy :data:`~ legacy.RingedPlanet` assembly helper.

        Returns
        -------
        :
            output : :data:`~ legacy.RingedPlanet`

        Caution
        --------
        It only works on single-planet system for improving performance
        """
        #Check if observer was provided
        if "Observer" in self.bodies:
            lamb=self.bodies["Observer"].lamb
            beta=self.bodies["Observer"].beta
            
        physics_defaults=deepcopy(LEGACY_PHYSICAL_PROPERTIES)
        physics_defaults.update(dict(limb_cs=self.bodies["Star"].limb_coeffs))
        physics_defaults.update(physics)
    
        #--CONSISTENCY--
        self._ringedplanet=dict(
            
            #Behavior
            behavior=dict(shadows=True),
            
            #Units
            CU=CanonicalUnits(UL=self.ul,UM=self.um),
    
            #Basic
            Rstar=self.bodies["Star"].radius,
            Rplanet=self.bodies["Planet"].radius,
    
            Rint=self.bodies["Ring"].fi,
            Rext=self.bodies["Ring"].fe,
            i=self.bodies["Ring"].i,
    
            a=self.bodies["Planet"].a,e=self.bodies["Planet"].e,
    
            #Orbit 
            Mstar=1,x=0,lambq=0,t0=0,kepler=False,
    
            #Observer
            eobs_ecl=np.array([lamb,beta]),
    
            #Sampling
            Np=self.bodies["Planet"].nspangles,
            Nr=self.bodies["Ring"].nspangles,
    
            Nb=0,Ns=30,
    
            #Physical properties
            physics=physics_defaults,
        )
        self.RP=RingedPlanet(**self._ringedplanet)
        return self.RP
    

    def update_StellarFlux(self):
        """
        Compute the Incident Stellar Flux for illuminated spangles that do not belong
        to stars, considering the effective area per spangle and applying the Flux Law

        Attributes
        ----------
        stellar_flux : `pd.Series`
            The computed incident stellar flux, stored in the Spangles Data object
            (``self.data.stellar_flux``).

        Note
        ---------
        - Only illuminated spangles not associated with stars are considered in the computation.
        - The computation follows the Flux Law (:math:`\propto d^{-2}`)
        
        .. math::
            F = \\frac{L}{4 \\pi d^2}
        """
        #Considered Spangles
        body_names = [name for name, body in self.bodies.items() if body.kind != 'Star']

        cond = self.data.name.isin(body_names)*self.data.illuminated

        #Creating stellar_flux attribute
        self.data['stellar_flux'] = np.zeros(self.data.shape[0])

        #Computing Incident Stellar Flux
        self.data.stellar_flux[cond] = abs(self.data.asp_luz[cond]*self.data.cos_luz[cond]/(4*np.pi*self.data.d_luz[cond]**2))


    def update_DiffuseReflection(self):
        """
        Compute the Diffuse Reflected Stellar Flux per Spangle that are both visible and illuminated,
        taking into account the illuminated side of the spangles and applying Lambert's Cosine Law.

        Returns
        ------------
        reflected_flux : `pd.Series`
            The computed reflected flux, stored in the Spangles Data object
            (``self.data.reflected_flux``).

        Note
        -------
        - Only visible and illuminated spangles reflect the incident stellar flux.
        - Diffuse reflection considers the illuminated side of the spangles, where the condition ``cos_obs * cos_luz > 0`` ensures the observer perceives the illuminated side.
        - The computation follows Lambert's Cosine Law (Lommel - Selliger Law is also available, see :any:`scatterer.Scatterer` for details).

        .. math::
            \\frac{F}{F_0} = \\sum_i \\frac{a_{s_i}\\cos \\Lambda_i}{4 \\pi d_i^2}A_{L_i}(\\Lambda_i)\\cos Z_i

        where:

        - :math:`F/F_0` is the reflected flux.
        - :math:`a_{s_i}` is the area of the spangle.
        - :math:`\\Lambda_i` is the angle between the incident light and the surface normal.
        - :math:`d_i` is the distance from the spangle to the light source.
        - :math:`A_{L_i}(\\Lambda_i)` is the Lambertian albedo of the spangle at angle :math:`\\Lambda_i`.
        - :math:`Z_i` is the angle between the observer's line of sight and the surface normal.
        """
        #Considered Spangles
        cond = self.data.illuminated*self.data.visible*(self.data.cos_obs*self.data.cos_luz > 0)*(~self.data.hidden)

        #Creating reflected_flux attribute
        self.data['reflected_flux'] = np.zeros(self.data.shape[0])

        #Compute incident stellar flux
        self.update_StellarFlux()
        
        #Update Lambertian Albedo Values for Integration Step
        self.update_albedos()

        #Computing Diffuse Reflected Light
        self.data.reflected_flux[cond] = self.data.stellar_flux[cond]*self.data.lambertian_albedo[cond]*self.data.cos_obs[cond]


    def update_Transit(self):
        """
        This method calculates the stellar flux drop caused by transiting spangles over star-kind
        bodies, accounting for limb-darkening effects and the spangle's optical properties.

        Returns
        ----------
        transit_flux : `pd.Series`
            The computed stellar flux drop, stored in the Spangles Data object
            (``self.data.transit_flux``).

        See Also
        ------------
        - The computation follows limb-darkening laws, as described in:
        `<https://pages.jh.edu/~dsing3/David_Sing/Limb_Darkening.html>`_.
        - Limb-darkening coefficients are sourced from:
        `<https://pages.jh.edu/~dsing3/LDfiles/LDCs.CoRot.Table1.txt>`_.

        Note
        ----------
        - Only spangles transiting over star-kind bodies are considered.
        - The flux drop is calculated as:

        .. math::
            \\frac{F}{F_0}= \\sum_i a_{s_i} \\cos Z_i \\beta_i(Z_i) \\frac{I(\\mu)}{I_0}

        where:

        - :math:`F/F_0` is the normalized stellar flux drop.
        - :math:`a_{s_i}` is the area of the spangle.
        - :math:`Z_i` is the angle between the observer's line of sight and the surface normal.
        - :math:`\\beta_i(Z_i)` is the attenuation factor for the spangle at angle :math:`Z_i` **[1]**
        .. math::
            \\beta_i(Z_i) = 1 - \\frac{\\tau}{2\\cos Z_i}e^{-\\frac{\\tau}{\\cos Z_i}}
        - :math:`I(\\mu)/I_0` is the intensity of the light at projected distance over stellar disk (see :any:`science.Science.limb_darkening` for details).
        
        **[1]**  French, R.G., Nicholson, P.D., 2000. Icarus 145, 502–523. doi:10. 1006/icar.2000.6357. 

        """
        #Creating transit_flux attribute
        self.data['transit_flux'] = np.zeros(self.data.shape[0])

        #Over all Stars
        for star in self.bodies:

            if self.bodies[star].kind == 'Star':

                #Considered Spangles
                cond = self.data.transit_over_obs.str.fullmatch(rf'^{star}[^&]*&$')*(~self.data.hidden)

                #Source-Spangle Projected Distance
                rhos = self.data.transit_over_obs[cond].str.extract(r':([^:]*)&', expand = False).astype(float)

                #Optical Parameters
                limb_coeffs = self.bodies[star].limb_coeffs
                norm_limb_coeff = self.bodies[star].norm_limb_darkening

                beta_values = 1 - np.exp(-self.data.tau_gray_optical[cond]/abs(self.data.cos_obs[cond]))

                star_scale = self.bodies[star].radius
                # limb_darkening = Util.limbDarkening(rhos, 1, limb_coeff, norm_limb_coeff)
                limb_darkening = Science.limb_darkening(rhos, cs = limb_coeffs , N = norm_limb_coeff)

                #Computing Stellar Flux Drop
                self.data.transit_flux[cond] += beta_values*abs(self.data.cos_obs[cond])*limb_darkening*(self.data.asp_obs[cond]/star_scale**2)
