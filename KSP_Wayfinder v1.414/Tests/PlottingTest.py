# -*- coding: utf-8 -*-
"""
Created on Wed 03.11.21
@author: v.fave
"""

"""for the import a dir above"""
import sys
sys.path.append("../WayfinderCore")

"""test the plotting functions"""
def run_plottingTest():

    from _Wayfinder import Wayfinder 

    plans = Wayfinder(planet_pack = "Vanilla")
    

    swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]]
    plans.load_df(datastore_name= "PlottingTest")  
    plans.find_best_plan(swing_by_bodies,t0_range=[0,1000])
    plans.plot_by_sequences(swing_by_bodies)
    plans.plot_DVvsT0(swing_by_bodies)

''' This is required to avoid issues with mutiprocessing '''
if __name__ == "__main__":
    run_plottingTest()