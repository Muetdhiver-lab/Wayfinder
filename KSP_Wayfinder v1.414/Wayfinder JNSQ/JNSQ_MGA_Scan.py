# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 12:38:45 2020

@author: v.fave
"""
import sys
sys.path.append("../WayfinderCore")

from _Wayfinder import Wayfinder 

plans = Wayfinder(planet_pack = "JNSQ")

setup = False

if setup :
    plans.add_batch(
            #swing_by_bodies = [[Kerbin],[Eve,Duna],[Kerbin,Eve,"*"],[Kerbin,Eve,"*"],[Moho]],
            swing_by_bodies = [["Kerbin"],["Kerbin"],["Kerbin"],["Dres"]],
            t0_min      = 0,
            t0_bin      = 200,
            n_t0_bins   = 5,
            auto_tof    = True,
            opt_level   = "debug",
            overwrite   = False)


else :
    swing_by_bodies =  [["Kerbin"],["Kerbin"],["Kerbin"],["Dres"]]
    plans.load_df()  


    #plansB.debugPrint()
    
    #plansB.edit_batch(swing_by_bodies, action = 'high')
    #plansB.optimize(n=5000)
    plans.optimize(n=20)
    plans.find_best_plan(swing_by_bodies)
    plans.plot_by_sequences(swing_by_bodies)

    #plansB.decode_solutions([["Kerbin"],["Eve"],["Moho"]])
