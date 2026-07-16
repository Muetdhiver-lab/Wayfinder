# -*- coding: utf-8 -*-
"""
Created on Wed Nov  18 12:38:45 2021

@author: v.fave

Dres being hard to get and unloved, lets make getting there fancy.
This example will generate several sequence in combinatoric fashion.
The sequence explored will be : KKKDr, KEKDr, KEDDr and KKDDr
T0 will be from 0 to 600 with 200 Kdy bins, Tof is automatic.
Injection set to circular, 100km. 
At the end of the jobs, the different sequence are plotted.
"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")
   
from _Wayfinder import Wayfinder 

DB_PATH = "wayfinder_examples.sqlite"
BATCH_NAME = "Ex3_DresAwareness"

def gen_vanillaDres_ex3():

    plans = Wayfinder(planet_pack = "Vanilla")
    
    swing_by_bodies = [["Kerbin"],["Kerbin","Eve"],["Kerbin","Duna"],["Dres"]]
    
    plans.add_batch_sqlite(
        db_path         = DB_PATH,
        batch_name      = BATCH_NAME,
        swing_by_bodies = swing_by_bodies,
        t0_min          = 0,
        t0_bin          = 200,
        n_t0_bins       = 3,          
        auto_tof        = True,
        opt_level       = "moderate",
        opt_injection   = "circular",
        injection_altitude = 100000)
    
def run_vanillaDres_ex3(): 


    plans = Wayfinder(planet_pack = "Vanilla")
    for i in range(4):
        plans.optimize_sqlite(DB_PATH, n=3, batch_name=BATCH_NAME)
    plans.find_best_known_plan_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        start_body="Kerbin",
        target_body="Dres",
        t0_range=[0,1000])
    plans.plot_by_sequences_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        start_body="Kerbin",
        target_body="Dres",
        t0_range=[0,1000])
    plans.plot_DVvsT0_sqlite(
        db_path=DB_PATH,
        batch_name=BATCH_NAME,
        start_body="Kerbin",
        target_body="Dres",
        t0_range=[0,1000])

''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    gen_vanillaDres_ex3()
    run_vanillaDres_ex3()
