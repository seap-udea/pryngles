import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pryngles import *
from pryngles import pixx
import time
#extension="pixx"
extension="cpixx"
import multiprocessing as mp
import os,sys,glob

def pool_handler(loc_num, params, func, multiplier):

    max_process = math.floor(mp.cpu_count()*multiplier)
    
    # Max allowed thread count on the server is 14
    if max_process>14: max_process=14
    
    if loc_num > max_process:
        run_num = 0
        while loc_num > 0:
            lbound = run_num*max_process
            if loc_num > max_process: ubound = (run_num+1)*max_process
            else: ubound = lbound + loc_num
            processes = [mp.Process(target=func, args=(params[i]))
                        for i in range(lbound, ubound)]
            for p in processes:
                time.sleep(0.1)
                p.start()
            for p in processes:
                time.sleep(0.1)
                p.join()
            loc_num = loc_num - max_process
            run_num+=1
    else:
        processes = [mp.Process(target=func, args=(params[i]))
                    for i in range(0, loc_num)]
        for p in processes:
            time.sleep(0.1)
            p.start()
        for p in processes:
            time.sleep(0.1)
            p.join()
        
def parametersweepGeom(ring_i_arr: np.ndarray,
                       fou_file_ring: str,
                       fou_file_planet: str,
                       orbit_i: float,
                       ring_l: float,
                       ring_ri: float,
                       ring_re: float,
                       tau_ring: float = 0.4,
                       Ns: int = 30,
                       Nb: int = 0,
                       Np: int = 10000,
                       Nr: int = 10000):
    
    # Announce start of function
    print (f"\n\n start run with: orbit i = {orbit_i}, ring l = {ring_l}") 
    
    # Generate save location if necessary
    if not os.path.isdir(f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{ring_i_arr[0]}"):      
        os.makedirs(f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{ring_i_arr[0]}")
    save_location = f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{ring_i_arr[0]}/"

    # Save log file since no printing comes out of multiprocessing
    sys.stdout = open(save_location + "logfile.out", "w")
    
    # Make sure all output is printed
    Verbose.VERBOSITY=VERB_ALL
    
    # Calculate starting position of observer and star
    gamma, beta_obs, lamb_obs, lamb_star = Util.calcStartingPosition(orbit_i,ring_i_arr[0],ring_l)
      
    # Initialise the system
    pixx_sys = System()
    s=pixx_sys.add(kind="Star",physics=dict(radius=Consts.rsun/pixx_sys.ul),optics=dict(limb_coeffs=[0.65]))
    p=pixx_sys.add(kind="Planet", primary=s,
                   orbit=dict(a=1, e=0.0),
                   physics=dict(radius=Consts.rsaturn/pixx_sys.ul),
                   optics=dict(nspangles=Np))
    r=pixx_sys.add(kind="Ring", primary=p,
                   physics=dict(fi=ring_ri, fe=ring_re, i=gamma),
                   optics=dict(nspangles=Nr))
    
    RP=pixx_sys.ensamble_system(extension=extension,
                                fname_planet=fou_file_planet,
                                fname_ring=fou_file_ring)
    
    lamb_initial = lamb_star
    lamb_final = lamb_initial + 360*Consts.deg
    lambs = np.linspace(lamb_initial,lamb_final,361)
        
    # Start series looping over the given ring inclinations
    for jj,r_i in enumerate(ring_i_arr):
        # Re-initialise the system for a different ring inclination
        if jj > 0:
            gamma, beta_obs, lamb_obs, lamb_star = Util.calcStartingPosition(orbit_i,r_i,ring_l)
            RP.i = gamma
            RP.changeObserver([lamb_obs,beta_obs])
            RP.Ns = Ns
            RP.Np = Np
            RP.Nr = Nr
            RP.Nb = Nb
            RP.updateProperties()
            lamb_initial = lamb_star
            lamb_final = lamb_initial + 360*Consts.deg
            lambs = np.linspace(lamb_initial,lamb_final,361)
            
        print (f"\n\n start run with: orbit i = {orbit_i}, ring l = {ring_l}, ring i = {r_i}") 
        
        # Generate save location if necessary
        if not os.path.isdir(f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{r_i}"):      
            os.makedirs(f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{r_i}")
        save_location = f"/home/allard/Data/Geom_series/Orbit_i_{orbit_i}/Ring_L_{ring_l}/Ring_i_{r_i}/"
        
        # Save log file since no printing comes out of multiprocessing
        sys.stdout = open(save_location + "logfile.out", "w")
        
        # Initialise the starting position
        RP.changeObserver([lamb_obs,beta_obs])
        RP.changeStellarPosition(lamb_initial)
        RP._updateGeometricalFactors()
        RP._updateIncomingStellarFlux()
        RP._updateObservedFacetAreas()
        
        # Save images showing the starting position of planet, ring and star
        ecl_fig,obs_fig,star_fig = RP.plotRingedPlanet(showstar=True,showfig=False)
        ecl_fig.savefig(save_location + f"fig_with_oi_{orbit_i}_rl_{ring_l}_ri_{r_i}_ecl.png", dpi=300)
        obs_fig.savefig(save_location + f"fig_with_oi_{orbit_i}_rl_{ring_l}_ri_{r_i}_obs.png", dpi=300)
        star_fig.savefig(save_location + f"fig_with_oi_{orbit_i}_rl_{ring_l}_ri_{r_i}_star.png", dpi=300)
        plt.close()
        
        # Make lists
        Stot  = []
        Sp    = []
        Sr    = []
        Ptot  = []
        Pp    = []
        Pr    = []
        alpha = []

        # Start the orbit
        for lamb in lambs:
            RP.changeStellarPosition(lamb)
            print("True anomaly: ", (lamb-lamb_initial)/Consts.deg)
            RP._updateGeometricalFactors()
            RP._updateIncomingStellarFlux()
            RP._updateObservedFacetAreas()
            RP.updateReflection(taur=tau_ring)
            
            # Save the relevant data
            Stot  += [RP.Stot]
            Sp    += [RP.Stotp]
            Sr    += [RP.Stotr]
            Ptot  += [RP.Ptot]
            Pp    += [RP.Ptotp]
            Pr    += [RP.Ptotr]
            alpha += [np.arccos(RP.alphaps)/Consts.deg]
            
        true_anomaly = list((lambs-lamb_initial)/Consts.deg)
        save_dict = {"lambda": true_anomaly, "alpha": alpha, "Stot": Stot,
                     "Sp": Sp, "Sr": Sr, "Ptot": Ptot, "Pp": Pp, "Pr": Pr}
        
        # Pickle the data, if file already exists it will be overwritten
        with open(save_location + f"data_with_oi_{orbit_i}_rl_{ring_l}_ri_{r_i}.pkl", "wb") as f:
            pickle.dump(save_dict, f)
            
def parameterSweep(fou_file_ring: str,
                   fou_file_planet: str,
                   orbit_i: float,
                   ring_i: float,
                   ring_l: float,
                   ring_ri: float,
                   ring_re: float,
                   name: str,
                   value,
                   tau_ring: float = 0.4,
                   Ns: int = 30,
                   Nb: int = 0,
                   Np: int = 10000,
                   Nr: int = 10000):
    
    # Announce start of function
    print (f"\n\n start run with: {name} = {value}") 
    
    # Generate save location if necessary
    if not os.path.isdir(f"/home/allard/Data/{name}_Series/{name}_{value}"):      
        os.makedirs(f"/home/allard/Data/{name}_Series/{name}_{value}")
    save_location = f"/home/allard/Data/{name}_Series/{name}_{value}/"

    # Save log file since no printing comes out of multiprocessing
    sys.stdout = open(save_location + "logfile.out", "w")
    
    # Make sure all output is printed
    Verbose.VERBOSITY=VERB_ALL
    
    # Calculate starting position of observer and star
    gamma, beta_obs, lamb_obs, lamb_star = Util.calcStartingPosition(orbit_i,ring_i,ring_l)
      
    # Initialise the system
    pixx_sys = System()
    s=pixx_sys.add(kind="Star",physics=dict(radius=Consts.rsun/pixx_sys.ul),optics=dict(limb_coeffs=[0.65]))
    p=pixx_sys.add(kind="Planet", primary=s,
                   orbit=dict(a=1, e=0.0),
                   physics=dict(radius=Consts.rsaturn/pixx_sys.ul),
                   optics=dict(nspangles=Np))
    r=pixx_sys.add(kind="Ring", primary=p,
                   physics=dict(fi=ring_ri, fe=ring_re, i=gamma),
                   optics=dict(nspangles=Nr))
    
    RP=pixx_sys.ensamble_system(extension=extension,
                                fname_planet=fou_file_planet,
                                fname_ring=fou_file_ring)
    
    lamb_initial = lamb_star
    lamb_final = lamb_initial + 360*Consts.deg
    lambs = np.linspace(lamb_initial,lamb_final,361)
        
    # Initialise the starting position
    RP.changeObserver([lamb_obs,beta_obs])
    RP.changeStellarPosition(lamb_initial)
    RP._updateGeometricalFactors()
    RP._updateIncomingStellarFlux()
    RP._updateObservedFacetAreas()

    # Save images showing the starting position of planet, ring and star
    ecl_fig,obs_fig,star_fig = RP.plotRingedPlanet(showstar=True,showfig=False)
    ecl_fig.savefig(save_location + \
                    f"fig_with_{name}_{value}_and_oi_{orbit_i}_rl_{ring_l}_ri_{ring_i}_rin_{ring_ri}_rout_{ring_re}_ecl.png", dpi=300)
    obs_fig.savefig(save_location + \
                    f"fig_with_{name}_{value}_and_oi_{orbit_i}_rl_{ring_l}_ri_{ring_i}_rin_{ring_ri}_rout_{ring_re}_obs.png", dpi=300)
    star_fig.savefig(save_location + \
                     f"fig_with_{name}_{value}_and_oi_{orbit_i}_rl_{ring_l}_ri_{ring_i}_rin_{ring_ri}_rout_{ring_re}_star.png", dpi=300)
    plt.close()

    # Make lists
    Stot  = []
    Sp    = []
    Sr    = []
    Ptot  = []
    Pp    = []
    Pr    = []
    alpha = []
    
    print("\n########################################")
    start_msg = f" Starting orbit simulation with: \n fou_file_planet = {RP.fname_planet} ,"+\
                f"\n fou_file_ring = {RP.fname_ring} ,"+\
                f"\n orbit inclination = {orbit_i} ,"+\
                f"\n ring roll = {ring_l} ,"+\
                f"\n ring inclination = {ring_i} ,"+\
                f"\n ring ri = {RP.fi} ,"+\
                f"\n ring re = {RP.fe} ,"+\
                f"\n observer i = {RP.eobs_ecl[1]/Consts.deg} ,"+\
                f"\n observer longitude = {RP.eobs_ecl[0]/Consts.deg} ,"+\
                f"\n ring inclination wrt ecl = {RP.i/Consts.deg} ,"+\
                f"\n number of ring spangles = {RP.Nrt} ,"+\
                f"\n number of planet spangles = {RP.Np} ,"+\
                f"\n ring opacity = {tau_ring} , "+\
                f"\n reference plane = {RP.reference_plane} "
    print(start_msg)
    print("########################################\n")
    
    # Start the orbit
    for lamb in lambs:
        RP.changeStellarPosition(lamb)
        print("True anomaly: ", (lamb-lamb_initial)/Consts.deg)
        RP._updateGeometricalFactors()
        RP._updateIncomingStellarFlux()
        RP._updateObservedFacetAreas()
        RP.updateReflection(taur=tau_ring)
        print("used ring opacity: ", RP.taur)
        
        # Save the relevant data
        Stot  += [RP.Stot]
        Sp    += [RP.Stotp]
        Sr    += [RP.Stotr]
        Ptot  += [RP.Ptot]
        Pp    += [RP.Ptotp]
        Pr    += [RP.Ptotr]
        alpha += [np.arccos(RP.alphaps)/Consts.deg]

    true_anomaly = list((lambs-lamb_initial)/Consts.deg)
    save_dict = {"lambda": true_anomaly, "alpha": alpha, "Stot": Stot,
                 "Sp": Sp, "Sr": Sr, "Ptot": Ptot, "Pp": Pp, "Pr": Pr}

    # Pickle the data, if file already exists it will be overwritten
    with open(save_location + f"data_with_{name}_{value}.pkl", "wb") as f:
        pickle.dump(save_dict, f)
        
def opticalThicknessTest(fou_file_num,orbit_i_arr,ring_i,ring_l):
    #Global parameters
    Ns=30
    Nb=0
    Np=500
    Nr=2000
    
    print ('\n\n start file: ', fou_file_num)  
    
    gamma, beta_obs, lamb_obs, lamb_star = Util.calcStartingPosition(orbit_i_arr[0],ring_i,ring_l)
    fou_ring = "./fou_files/Ring/fou_ring_"+fou_file_num+"_0_8.dat"
    
    sys_test = System()
    s=sys_test.add(kind="Star",physics=dict(radius=Consts.rsun/sys_test.ul),optics=dict(limb_coeffs=[0.65]))
    p=sys_test.add(kind="Planet", primary=s,
                   orbit=dict(a=3, e=0.0),
                   physics=dict(radius=Consts.rsaturn/sys_test.ul),
                   optics=dict(nspangles=Np))
    r=sys_test.add(kind="Ring", primary=p,
                   physics=dict(fi=7, fe=9.25, i=gamma),
                   optics=dict(nspangles=Nr))
    RP_test=sys_test.ensamble_system(extension=extension,
                                     fname_planet=Misc.get_data("fou_lambert.dat"),
                                     fname_ring=fou_ring)
    Rrs = np.zeros((len(orbit_i_arr),3))
    Prs = np.zeros(len(orbit_i_arr))
    illum_angle = np.zeros(len(orbit_i_arr))
    view_angle = np.zeros(len(orbit_i_arr))

    for jj,o_i in enumerate(orbit_i_arr):
        if jj > 0:
            gamma, beta_obs, lamb_obs, lamb_star = Util.calcStartingPosition(o_i,ring_i,ring_l)
            RP_test.i = gamma
            RP_test.changeObserver([lamb_obs,beta_obs])
            RP_test.Ns = Ns
            RP_test.Np = Np
            RP_test.Nr = Nr
            RP_test.Nb = Nb
            RP_test.updateProperties()
        
        RP_test.changeObserver([lamb_obs,beta_obs])
        RP_test.changeStellarPosition(lamb_star)
        RP_test.updateOpticalFactors()
        RP_test.updateReflection()
        
        if (np.inner(RP_test.nstar_equ,RP_test.nr_equ) < 0) ^ (np.inner(RP_test.nobs_equ,RP_test.nr_equ) < 0):
            Rrs[jj] = RP_test.Stotr
            Prs[jj] = RP_test.Ptotr
            
        illum_angle[jj] = RP_test.etars[0]
        view_angle[jj] = RP_test.zetars[0]
    
    save_dict = {"orbit_i": orbit_i_arr, "illum": illum_angle, "Flux": Rrs, "Degree": Prs}
    with open(f"./data/optical_thickness/opt_thickness_{fou_file_num}.pkl", "wb") as f:
        pickle.dump(save_dict, f)