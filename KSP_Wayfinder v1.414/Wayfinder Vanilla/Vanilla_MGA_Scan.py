# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 12:38:45 2020

@author: v.fave
"""
import sys
sys.path.append("../WayfinderCore")

def run_vanilla_mga_scan():

    from _Wayfinder import Wayfinder 

    plans = Wayfinder(planet_pack = "Vanilla")
    
    setup = False
    
    if setup :
        plans.add_batch(
            #swing_by_bodies = [[Kerbin],[Eve,Duna],[Kerbin,Eve,"*"],[Kerbin,Eve,"*"],[Moho]],
            swing_by_bodies = [["Kerbin"],["Eve"],["Kerbin"],["Kerbin"],["Jool"],["Eeloo"]],
            #swing_by_bodies = [["Kerbin"],["Eve","Kerbin"],["Kerbin","Duna"],["Dres"]],
            t0_min          = 0,
            t0_bin          = 200,
            n_t0_bins       = 20,          
            auto_tof        = True,
            opt_level       = "ultra+",
            opt_injection   = "circular",
            injection_altitude = 100000)



    else :
        swing_by_bodies = [["Kerbin"],["Eve"],["Moho"]]
        #swing_by_bodies = [["Kerbin"],["Eve"],["Kerbin"],["Kerbin"],["Jool"],["Eeloo"]]
        #swing_by_bodies = [["Kerbin"],["Eve","Kerbin"],["Kerbin","Duna"],["Dres"]]
        #swing_by_bodies = "KKEMo"
        plans.load_df()  
        #print(plans.generateSequences(swing_by_bodies))
        plans.auto_tof(["Kerbin","Kerbin","Duna","Dres"],debug=True)
        #plans.recalc_results()
        plans.audit_results()
        #plans.edit_batch(swing_by_bodies, action = "inj:vinf")
        #plans.edit_batch(swing_by_bodies, action = "lvl:high")

        plans.find_best_plan(swing_by_bodies,t0_range=[000,2000])
        #plans.plot_by_sequences(swing_by_bodies,t0_range=[0,2000])
        #plans.plot_DVvsT0(swing_by_bodies,t0_range=[0,2000])
        #plans.decode_solutions(swing_by_bodies)

''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    run_vanilla_mga_scan()