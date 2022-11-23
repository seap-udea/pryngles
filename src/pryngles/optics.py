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
#!/usr/bin/env python
# coding: utf-8

# # Pryngles module: optics
# 
# Template of a module

# ## External modules

#@external
from pryngles import *
#@end:external
import pryngles.pixx as pixx
import time
pixx.rdfous_planet.__doc__

Scatterer_doc="""This is the basic class of a scatterer
"""

# ## Class: Scatterer

#@class
class Scatterer(PrynglesCommon):
    def __init__(self,
                 fname_planet: str = None,
                 fname_ring: str = None):
        """
        """
        self.save_values=[]
        
        self.data = None
        self.SPANGLER_SCATTERER_COLUMNS_UPDATE = ["name","x_obs","y_obs","z_obs","asp",
                                                  #Angles
                                                  "cos_luz","cos_obs","azim_obs_luz",
                                                  #State
                                                  "visible","shadow","indirect","emit","illuminated","transmit",
                                                  "hidden_by_luz","transit_over_luz","hidden_by_obs","transit_over_obs",
                                                  #Transit
                                                  "transit","occult",
                                                  #Physical properties
                                                  "albedo_gray_normal", "tau_gray_optical"]
        
        # Read scatter data,
        # !!! The filenames need to be less than 100 characters !!!
        self.read_data(fname_planet, fname_ring)
        
        # Set number of stokes elements that need to be calculated
        if fname_planet is None and fname_ring is None:
            self.nmatp = 1
            self.nmatr = 1
            self.nmat = 1
        elif fname_planet is None:
            self.nmatp = 1
            self.nmat = self.nmatr
        elif fname_ring is None:
            self.nmatr = 1
            self.nmat = self.nmatp
        else:
            self.nmat = np.max([self.nmatp,self.nmatr])
        
        # Output
        self.Stotp = np.zeros(self.nmatp)
        self.Ptotp = 0 
        self.Stotr = np.zeros(self.nmatr)
        self.Ptotr = 0 
        self.Stot = np.zeros(self.nmat)
        self.Ptot = 0
        
        # Add columns to dataframe
        if self.nmat == 4:
            if self.nmatr == 4:
                self.STOKES_VECTOR_RING = ["F","Q","U","V","P"]
            else:
                self.STOKES_VECTOR_RING = ["F","Q","U","P"]
            if self.nmatp == 4: 
                self.STOKES_VECTOR_PLANET = ["F","Q","U","V","P"]
            else:
                self.STOKES_VECTOR_PLANET = ["F","Q","U","P"]
            self.SPANGLER_SCATTERER_COLUMNS = ["beta_loc","F","Q","U","V","P"]
        elif self.nmat == 3:
            self.STOKES_VECTOR_RING = ["F","Q","U","P"]
            self.STOKES_VECTOR_PLANET = ["F","Q","U","P"]
            self.SPANGLER_SCATTERER_COLUMNS = ["beta_loc","F","Q","U","P"]
        elif self.nmat == 1:
            self.STOKES_VECTOR_RING = ["F"]
            self.STOKES_VECTOR_PLANET = ["F"]
            self.SPANGLER_SCATTERER_COLUMNS = ["beta_loc","F"]
        
    def read_data(self, 
                 fname_planet: str,
                 fname_ring: str):
        """
        Reads-in the fourier coefficients from the specified files
        Reading is by a FORTRAN function
        """
        if fname_planet is not None:
            i = j = 0
            with open(fname_planet) as file:
                for line in file: 
                    if line.rstrip()[0] != "#":
                        if i == 0:
                            nmatp = int(line.rstrip())
                        elif i ==1:
                            nmugsp = int(line.rstrip())
                        else:
                            j += 1
                        i += 1

            nfoup = int( (j-nmugsp)/(nmugsp**2) )
            
            self.nfoup = nfoup
            self.nmatp = nmatp
            self.nmugsp = nmugsp
            
            # Reflected light
            self.xmup,self.rfoup = pixx.rdfous_planet(fname_planet,nfoup,nmatp,nmugsp)
                        
        if fname_ring is not None:
            i = j = 0
            with open(fname_ring) as file:
                for line in file: 
                    if line.rstrip()[0] != "#":
                        if i == 0:
                            nmatr = int(line.rstrip())
                        elif i ==1:
                            nmugsr = int(line.rstrip())
                        else:
                            j += 1
                        i += 1

            nfour = int( (j-nmugsr)/(nmugsr**2) )
            
            self.nfour = nfour
            self.nmatr = nmatr
            self.nmugsr = nmugsr
            
            # Reflected light
            self.xmur,self.rfour = pixx.rdfous_ring(fname_ring,False,nfour,nmatr,nmugsr)
            
            # Transmitted light
            self.xmur,self.tfour = pixx.rdfous_ring(fname_ring,True,nfour,nmatr,nmugsr)
        
    def update_data(self,system):
        """
        Copies select data from the system object to a dataframe
        Adds the columns containing the beta angle, Stokes vector and degree of polarization
        """
        self.sys = system
        cond = system.sg.data.name != "Star"
        if self.data is None:
            self.data = deepcopy(system.sg.data.loc[cond,self.SPANGLER_SCATTERER_COLUMNS_UPDATE])
        else:
            self.data[self.SPANGLER_SCATTERER_COLUMNS_UPDATE] =                 deepcopy(system.sg.data.loc[cond,self.SPANGLER_SCATTERER_COLUMNS_UPDATE])
        self.data[self.SPANGLER_SCATTERER_COLUMNS] = np.zeros((cond.sum(),len(self.SPANGLER_SCATTERER_COLUMNS)))
        
        # Reset output
        self.Stotp = np.zeros(self.nmatp)
        self.Ptotp = 0 
        self.Stotr = np.zeros(self.nmatr)
        self.Ptotr = 0 
        self.Stot = np.zeros(self.nmat)
        self.Ptot = 0
        
    def compute_angles(self):
        """
        Function that computes:
            - The phase angle
            - The azimuthal difference angle for every spangle
            - The beta angle for every spangle
        
        The beta angle rotates the local scattering plane to the planetary scattering plane
        """
        for name,body in self.sys.bodies.items():
            if body.kind == "Planet":
                center = body.center_ecl
                        
        azim,incli = Science.spherical(self.sys.n_obs)[1:]
        Rx = self.rotation_matrix_x(np.pi/2-incli)
        Rz = self.rotation_matrix_z(np.pi/2-azim)
        
        # Calculate normal vector and distance to star
        luz_ecl,self.d_luz = spy.unorm(center-self.sys.center_root)
        self.luz_obs = np.matmul(Rx, np.matmul(Rz, luz_ecl))
        
        self.phase_angle = np.dot(self.luz_obs,np.array([0,0,1]))
        
        for name,body in self.sys.bodies.items():
            if body.kind == "Star":
                verbose(VERB_SIMPLE,f"Body is a star... skipping")
                continue
                
            elif body.kind == "Planet":
                cond = self.data.name == body.kind
                etaps = self.data.loc[cond,"cos_luz"]
                zetaps =  self.data.loc[cond,"cos_obs"]
                
                # Azimuthal angle difference 
                t1 = self.phase_angle - zetaps*etaps
                t2 = np.sin(np.arccos(etaps))*np.sin(np.arccos(zetaps))
                t3 = t1/t2
                t3[t3 > 1] = 1.0
                t3[t3 < -1] = -1.0
                phidiffps = np.pi - np.arccos(t3)
                phidiffps[abs(t2) <= 1e-9] = 0.0 
                phidiffps[self.data.loc[cond, "y_obs"] < 0] *= -1
                self.data.loc[cond,"azim_obs_luz"] = phidiffps
                
                # Calculate beta angle
                x = self.data.loc[cond,"x_obs"]
                y = self.data.loc[cond,"y_obs"]
                if self.luz_obs[0] >= 0:
                    betaps = np.arctan(y/x)
                    betaps[x*y < 0] += np.pi
                else:
                    betaps = -np.arctan(y/x)
                    betaps[x*y >= 0] += np.pi
                self.data.loc[cond,"beta_loc"] = betaps
                
            elif body.kind == "Ring":
                cond = self.data.name == body.kind
                etars = self.data.loc[cond,"cos_luz"]
                zetars =  self.data.loc[cond,"cos_obs"]
                
                # Azimuthal angle difference 
                t1 = self.phase_angle - zetars*etars
                t2 = np.sin(np.arccos(etars))*np.sin(np.arccos(zetars))
                t3 = t1/t2
                t3[t3 > 1] = 1.0
                t3[t3 < -1] = -1.0
                phidiffrs = np.pi - np.arccos(t3)
                phidiffrs[abs(t2) <= 1e-9] = 0.0 
                self.data.loc[cond,"azim_obs_luz"] = phidiffrs
                
                # Beta calculation for rings
                t1 = etars - self.phase_angle*zetars
                t2 = np.sin(np.arccos(self.phase_angle))*np.sin(np.arccos(zetars))
                t3 = t1/t2
                t3[t3 > 1] = 1.0
                t3[t3 < -1] = -1.0
        
                if self.luz_obs[0] >= 0:
                    betars = np.arccos(t3)
                    betars[abs(t2) <= 1e-9] = 0.0 
                else:
                    betars = np.pi - np.arccos(t3)
                    betars[abs(t2) <= 1e-9] = 0.0 
                self.data.loc[cond,"beta_loc"] = betars
                
            else: 
                continue
    
    def scattering(self,
                   system,
                   normalize: bool = False):
        """
        Function that first updates the local dataframe with new geometry data
        Then computes the angles needed for the scattering code to work
        
        Checks which spangles are active and for which the scattered stokes vector needs to calculated
        Passes all the necessary angles and the fourier-coefficients file to a FORTRAN code 
        In the FORTRAN code the fourier coefficients are interpolated to the required angles
        and then summed to calculate the stokes vector.
        """
        # Update positional data
        self.update_data(system)
        
        # Calculate necessary angles
        self.compute_angles()

        angle_eps = 1e-3 # Value above which Cos(angle) needs to be
        rings_present = False
        
        #Planet
        for name,body in self.sys.bodies.items():
            if body.kind == "Star":
                verbose(VERB_SIMPLE,f"Body is a star... skipping")
                continue       
            elif body.kind == "Planet":
                condp = self.data.name == body.kind
                normp = 4*np.pi*body.radius**2
                if body.childs:
                    for bhash,child in body.childs.items():
                        if child.kind == "Ring":
                            rings_present = True
            elif body.kind == "Ring":
                condr = self.data.name == body.kind
                normr = 4*np.pi*(body.re**2 - body.ri**2)
        
        for name,body in self.sys.bodies.items():
            if body.kind == "Planet":
                # Spangles that are visible, illuminated and not in transit
                condv = (self.data.loc[condp,"visible"])&(self.data.loc[condp,"illuminated"])&                        (~self.data.loc[condp,"transit"])
                
                # Check if there are rings present
                if rings_present:
                    # Facets that are illuminated through the rings
                    condspr = (self.data.hidden_by_luz.apply(lambda x:"Ring" in x))&(self.data.visible)&                              (~self.data.loc[condp,"transit"])

                    # Facets that are visible but the line of sight is blocked by the rings
                    condspo = (self.data.hidden_by_obs.apply(lambda x:"Ring" in x))&(self.data.illuminated)&                              (~self.data.loc[condp,"transit"])
                    
                    cond = (condv) | (condspr) | (condspo)
                else:
                    cond = np.logical_and(condv,condp)
                
                # Only proceed if there are active spangles
                if cond.sum() > 0:
                    
                    self.data.loc[cond,self.STOKES_VECTOR_PLANET] = pixx.reflection(cond.sum(),
                                                                             self.data.loc[cond,"azim_obs_luz"], 
                                                                             self.data.loc[cond,"beta_loc"],
                                                                             abs(self.data.loc[cond,"cos_luz"]), 
                                                                             abs(self.data.loc[cond,"cos_obs"]),
                                                                             self.nmugsp,self.nmatp,self.nfoup,
                                                                             self.xmup,self.rfoup,
                                                                             self.data.loc[cond,"asp"]/normp)
                    
                    self.save_values+=[
                        dict(
                            obj="planet",
                            args=[
                                cond.sum(),
                                self.data.loc[cond,"azim_obs_luz"].values, 
                                self.data.loc[cond,"beta_loc"].values,
                                abs(self.data.loc[cond,"cos_luz"]).values, 
                                abs(self.data.loc[cond,"cos_obs"]).values,
                                self.nmugsp,self.nmatp,self.nfoup,
                                self.xmup,self.rfoup,
                                self.data.loc[cond,"asp"].values/normp,
                            ],
                            stokes=self.data.loc[cond,self.STOKES_VECTOR_PLANET].values,
                        )
                    ]

                    # Check if the rings are seen edge-on and/or illuminated edge-on
                    if rings_present:
                        vcheck = abs(self.data.loc[condr,"cos_obs"].iloc[0]) > angle_eps # seen
                        icheck = np.mean(abs(self.data.loc[condr,"cos_luz"])) > angle_eps # illuminated
                        vsum = condspo.sum()
                        isum = condspr.sum()
                        if vcheck and icheck and (vsum > 0) and (isum > 0):
                            self.data.loc[condspr,self.STOKES_VECTOR_PLANET[:-1]] *=                                         np.exp(-abs(np.mean(self.data.loc[condr,"tau_gray_optical"]/
                                                            self.data.loc[condr,"cos_luz"])))
                            self.data.loc[condspo,self.STOKES_VECTOR_PLANET[:-1]] *=                                         np.exp(-self.data.loc[condr,"tau_gray_optical"].iloc[0]/
                                                      abs(self.data.loc[condr,"cos_obs"].iloc[0]))
                        elif vcheck and (vsum > 0):
                            self.data.loc[condspo,self.STOKES_VECTOR_PLANET[:-1]] *=                                         np.exp(-self.data.loc[condr,"tau_gray_optical"].iloc[0]/
                                                      abs(self.data.loc[condr,"cos_obs"].iloc[0]))
                        elif icheck and (isum > 0):
                            self.data.loc[condspr,self.STOKES_VECTOR_PLANET[:-1]] *=                                         np.exp(-abs(np.mean(self.data.loc[condr,"tau_gray_optical"]/
                                                            self.data.loc[condr,"cos_luz"])))

                    Stotp = np.sum(self.data.loc[cond,self.STOKES_VECTOR_PLANET[:-1]],axis=0)
                    if Stotp[0] < 1e-6:
                        Ptotp = 0
                    elif abs(Stotp[2]) < 1e-6:
                        Ptotp = -Stotp[1]/Stotp[0]
                    else:
                        if self.nmatp == 4:
                            Ptotp = np.sqrt(Stotp[1]**2 + Stotp[2]**2 + Stotp[3]**2)/Stotp[0]
                        else:
                            Ptotp = np.sqrt(Stotp[1]**2 + Stotp[2]**2)/Stotp[0]
                    
                    # The normalized stokes vectors are given units
                    if not normalize:
                        self.data.loc[cond,self.STOKES_VECTOR_PLANET[:-1]] /= (4*np.pi*self.d_luz**2)*1e6 #ppm
                        self.data.loc[cond,self.STOKES_VECTOR_PLANET[:-1]] *= normp
                    
                        # Set integrated stokes vector
                        self.Stotp = Stotp/(4*np.pi*self.d_luz**2)*1e6*normp 
                    else:
                        self.data.loc[cond,self.STOKES_VECTOR_PLANET[:-1]] *= 4
                        self.Stotp = Stotp*4
                        
                    self.Ptotp = Ptotp  
                    
            elif body.kind == "Ring":
                # Spangles that are visible and illuminated
                condv = (self.data.loc[condr,"visible"])&(self.data.loc[condr,"illuminated"])

                # Check if there is transmission through the ring
                transmission = False 
                if (np.mean(self.data.loc[condr,"cos_luz"]) < 0) ^ (self.data.loc[condr,"cos_obs"].iloc[0] < 0):
                    transmission = True

                # Make sure the dimensions are correct
                cond = np.logical_and(condv,condr)

                if cond.sum() > 0:                    
                    if transmission:
                        self.data.loc[cond,self.STOKES_VECTOR_RING] = pixx.reflection(cond.sum(),
                                                                               self.data.loc[cond,"azim_obs_luz"], 
                                                                               self.data.loc[cond,"beta_loc"],
                                                                               abs(self.data.loc[cond,"cos_luz"]), 
                                                                               abs(self.data.loc[cond,"cos_obs"]),
                                                                               self.nmugsr,self.nmatr,self.nfour,
                                                                               self.xmur,self.tfour,
                                                                               self.data.loc[cond,"asp"]/normr)
                        self.save_values+=[
                            dict(
                                obj="planet",
                                args=[
                                    cond.sum(),
                                    self.data.loc[cond,"azim_obs_luz"].values, 
                                    self.data.loc[cond,"beta_loc"].values,
                                    abs(self.data.loc[cond,"cos_luz"]).values, 
                                    abs(self.data.loc[cond,"cos_obs"]).values,
                                    self.nmugsr,self.nmatr,self.nfour,
                                    self.xmur,self.tfour,
                                    self.data.loc[cond,"asp"].values/normr
                                ],
                                stokes=self.data.loc[cond,self.STOKES_VECTOR_PLANET].values,
                            )
                        ]
                        
                    else:
                        self.data.loc[cond,self.STOKES_VECTOR_RING] = pixx.reflection(cond.sum(),
                                                                               self.data.loc[cond,"azim_obs_luz"], 
                                                                               self.data.loc[cond,"beta_loc"],
                                                                               abs(self.data.loc[cond,"cos_luz"]), 
                                                                               abs(self.data.loc[cond,"cos_obs"]),
                                                                               self.nmugsr,self.nmatr,self.nfour,
                                                                               self.xmur,self.rfour,
                                                                               self.data.loc[cond,"asp"]/normr)
                        
                        self.save_values+=[
                            dict(
                                obj="planet",
                                args=[
                                       cond.sum(),
                                       self.data.loc[cond,"azim_obs_luz"].values, 
                                       self.data.loc[cond,"beta_loc"].values,
                                       abs(self.data.loc[cond,"cos_luz"]).values, 
                                       abs(self.data.loc[cond,"cos_obs"]).values,
                                       self.nmugsr,self.nmatr,self.nfour,
                                       self.xmur,self.rfour,
                                       self.data.loc[cond,"asp"].values/normr
                                ],
                                stokes=self.data.loc[cond,self.STOKES_VECTOR_PLANET].values,
                            )
                        ]

                    Stotr = np.sum(self.data.loc[cond,self.STOKES_VECTOR_RING[:-1]],axis=0)
                    if Stotr[0] < 1e-6:
                        Ptotr = 0
                    elif abs(Stotr[2]) < 1e-6:
                        Ptotr = -Stotr[1]/Stotr[0]
                    else:
                        if self.nmatr == 4:
                            Ptotr = np.sqrt(Stotr[1]**2 + Stotr[2]**2 + Stotr[3]**2)/Stotr[0]
                        else:
                            Ptotr = np.sqrt(Stotr[1]**2 + Stotr[2]**2)/Stotr[0]
                    # The normalized stokes vectors are given units
                    if not normalize:
                        self.data.loc[cond,self.STOKES_VECTOR_RING[:-1]] /= (4*np.pi*self.d_luz**2)*1e6 #ppm
                        self.data.loc[cond,self.STOKES_VECTOR_RING[:-1]] *= normr
                        
                        # Set integrated stokes vector
                        self.Stotr = Stotr/(4*np.pi*self.d_luz**2)*1e6*normr 
                    else:
                        self.data.loc[cond,self.STOKES_VECTOR_RING[:-1]] *= 4
                        self.Stotr = Stotr*4
                        
                    self.Ptotr = Ptotr
            else:
                continue
                
        # Calculate total flux and total degree of polarization
        if self.nmat > 1:            
            Stot = np.sum(self.data.loc[:,self.SPANGLER_SCATTERER_COLUMNS[1:-1]],axis=0)
            if abs(Stot[2]) < 1e-18:
                Ptot = -Stot[1]/Stot[0]
            else:
                if self.nmatr == 4:
                    Ptot = np.sqrt(Stot[1]**2 + Stot[2]**2 + Stot[3]**2)/Stot[0]
                else:
                    Ptot = np.sqrt(Stot[1]**2 + Stot[2]**2)/Stot[0]
            self.Stot = Stot
            self.Ptot = Ptot
        else:
            self.Stot = self.Stotp + self.Stotr
            self.Ptot = 0
            
        output_dict = {"Stot": self.Stot, "Ptot": self.Ptot, "Stotp": self.Stotp,
                       "Ptotp": self.Ptotp, "Stotr": self.Stotr, "Ptotr": self.Ptotr}
        return output_dict

    def rotation_matrix_x(self,angle):
        """
        Rotation matrix for a rotation around the x-axis
        """
        Rm = np.array([[1,0,0],[0,np.cos(angle), np.sin(angle)],[0,-np.sin(angle),np.cos(angle)]])
        return Rm
    
    def rotation_matrix_z(self,angle):
        """
        Rotation matrix for a rotation around the z-axis
        """
        Rm = np.array([[np.cos(angle), -np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
        return Rm
        
    def orbital_rotation_rate(self):
        """
        Calculate the orbital rotation rate of a circular orbit
        Used in converting the integration time to degrees of true anomaly
        """
        for name,body in self.sys.bodies.items():
            if name == "Star":
                ms = body.m
            elif name == "Planet":
                a = body.a
                mp = body.m
                
        n = 1/np.sqrt(a**3/(self.sys.G*(ms+mp)))
        return n
    
    def lambertian_test(self,alpha):
        """
        Simple, analytical model for the normalized reflected light coming of a lambertian planet
        """
        F = 2*(np.sin(alpha)+(np.pi-alpha)*np.cos(alpha))/3
        return F
        
Scatterer.__doc__=Scatterer_doc
#@end:class


# ### The end

