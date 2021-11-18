# -*- coding: utf-8 -*-
"""
Created on Wed Nov  18 12:38:45 2021

@author: v.fave
"""

"""for the import a dir above"""
import sys
import time
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

def gen_vanilla_mga_scan_test():

    plans = Wayfinder(planet_pack = "JNSQ")
    
    swing_by_bodies = [["Kerbin"],["Kerbin"],["Duna"],["Edna"]]
    
    plans.add_batch(
        datastore_name  = 'Wayfinder_JNSQ_Benchmarking',
        overwrite       = True,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 1800,
        t0_bin          = 200,
        n_t0_bins       = 1,          
        auto_tof        = True,
        opt_level       = "debug",
        opt_injection   = "circular",
        injection_altitude = 100000,
        lambert_max_revs   = 0)
    
def run_vanilla_mga_scan_test(): 


    plans = Wayfinder(planet_pack = "JNSQ")
    plans.load_df(datastore_name = 'Wayfinder_JNSQ_Benchmarking')  
    swing_by_bodies = [["Kerbin"],["Kerbin"],["Duna"],["Edna"]]
    plans.optimize(n=1)
    plans.find_best_plan(swing_by_bodies)


''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_vanilla_mga_scan_test()
    starttime = time.time()
    run_vanilla_mga_scan_test()
    elapsed_time_fl = (time.time() - starttime)
    print("time to run : "+str(elapsed_time_fl))