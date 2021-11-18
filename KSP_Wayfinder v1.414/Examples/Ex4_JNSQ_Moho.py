# -*- coding: utf-8 -*-
"""
Created on Wed Nov  18 12:38:45 2021

@author: v.fave

A simple example of JNSQ Job for moho with either KEMo or KEEMo using the "*" wildcard to an optional assist, 
Two bins in T0 and tof. The edit method is used to alter the jobs optimization level and orbital injection to v_inf. 
(that is : the optimizer will try to minimize injection speed at Moho SOI, not the actual injection burn cost)
"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

def gen_JNSQ_Moho_ex4():

    plans = Wayfinder(planet_pack = "JNSQ")
    
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]]
    
    plans.add_batch(
        datastore_name  = 'Ex4_JNSQ_Moho',
        overwrite       = True,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 0,
        t0_bin          = 200,
        n_t0_bins       = 2,          
        tof_min         = 300,
        tof_bin         = 150,
        n_tof_bins      = 2,   
        opt_level       = "high",
        opt_injection   = "circular",
        injection_altitude = 100000)
    
def run_JNSQ_Moho_ex4(): 


    plans = Wayfinder(planet_pack = "JNSQ")
    plans.load_df(datastore_name = 'Ex4_JNSQ_Moho')  
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]]
    plans.edit_batch(swing_by_bodies, action = "inj:vinf")
    plans.edit_batch(swing_by_bodies, action = "lvl:moderate")
    plans.optimize(n=8)
    plans.find_best_plan(swing_by_bodies,t0_range=[0,1000])
    plans.plot_by_sequences(swing_by_bodies,t0_range=[0,1000])

''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_JNSQ_Moho_ex4()
    run_JNSQ_Moho_ex4()