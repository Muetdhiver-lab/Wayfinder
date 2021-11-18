# -*- coding: utf-8 -*-
"""
Created on Wed Nov  18 12:38:45 2021

@author: v.fave

A simple example of KEEMo grav assist search in the first 1000 days
"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

def gen_vanilla_mga_ex1():

    plans = Wayfinder(planet_pack = "Vanilla")
    
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve"],["Moho"]]
    
    plans.add_batch(
        datastore_name  = 'Ex1_KEEMo_MGA',
        overwrite       = True,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 0,
        t0_bin          = 1000,
        n_t0_bins       = 1,          
        auto_tof        = True,
        opt_level       = "high",
        opt_injection   = "circular",
        injection_altitude = 100000)
    
def run_vanilla_mga_ex1(): 


    plans = Wayfinder(planet_pack = "Vanilla",datastore_name = 'Ex1_KEEMo_MGA')
    plans.load_df(datastore_name = 'Ex1_KEEMo_MGA')  
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve"],["Moho"]]
    plans.optimize(n=1)
    plans.find_best_plan(swing_by_bodies,t0_range=[0,1000])


''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_vanilla_mga_ex1()
    run_vanilla_mga_ex1()