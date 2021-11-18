# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 12:38:45 2020

@author: v.fave
"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

def gen_JNSQ_mga_scan_test():

    plans = Wayfinder(planet_pack = "JNSQ")
    
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]]
    
    plans.add_batch(
        datastore_name  = 'JNSQ_MGA_Scan_test',
        overwrite       = True,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 0,
        t0_bin          = 200,
        n_t0_bins       = 5,          
        auto_tof        = True,
        opt_level       = "high",
        opt_injection   = "circular",
        injection_altitude = 100000)
    
def run_JNSQ_mga_scan_test(): 


    plans = Wayfinder(planet_pack = "JNSQ",datastore_name = 'JNSQ_MGA_Scan_test')
    plans.load_df(datastore_name = 'JNSQ_MGA_Scan_test')  
    swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]]
    plans.edit_batch(swing_by_bodies, action = "lvl:debug")
    plans.optimize(n=10)
    plans.find_best_plan(swing_by_bodies,t0_range=[0,1000])


''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_JNSQ_mga_scan_test()
    run_JNSQ_mga_scan_test()