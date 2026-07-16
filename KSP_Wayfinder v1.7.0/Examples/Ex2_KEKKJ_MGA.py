# -*- coding: utf-8 -*-
"""
Created on Wed Nov  18 12:38:45 2021

@author: v.fave

A simple example of KEKKJ grav assist search with a launch window from Kdy 600 to 800
and a Tof of 3000 to 4000 Kdy. Injection in elliptical orbit with pe set at 1000 km

"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

DB_PATH = "wayfinder_examples.sqlite"
BATCH_NAME = "Ex2_KEKKJ_MGA"

def gen_vanilla_mga_ex2():

    plans = Wayfinder(planet_pack = "Vanilla")
    
    swing_by_bodies = [["Kerbin"],["Eve"],["Kerbin"],["Kerbin"],["Jool"]]
    
    plans.add_batch_sqlite(
        db_path         = DB_PATH,
        batch_name      = BATCH_NAME,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 600,
        t0_bin          = 200,
        n_t0_bins       = 1,          
        tof_min          = 2500,
        tof_bin          = 1000,
        n_tof_bins       = 1,        
        opt_level       = "ultra",
        opt_injection   = "elliptical",
        injection_altitude = 1000000)
    
def run_vanilla_mga_ex2(): 


    plans = Wayfinder(planet_pack = "Vanilla",datastore_name = BATCH_NAME)
    swing_by_bodies = [["Kerbin"],["Eve"],["Kerbin"],["Kerbin"],["Jool"]]
    plans.optimize_sqlite(DB_PATH, n=1, batch_name=BATCH_NAME)
    plans.find_best_known_plan_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        sequence_short_name="KEKKJ",
        t0_range=[0,1000])


''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_vanilla_mga_ex2()
    run_vanilla_mga_ex2()
