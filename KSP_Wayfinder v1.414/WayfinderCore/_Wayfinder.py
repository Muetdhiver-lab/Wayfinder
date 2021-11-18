# Launch window scanner - coarse porkchop.
"""
Created on Wed May 29 08:02:50 2019

@author: v.fave
"""
import os
import copy
import pandas as pd
import pygmo as pg
import numpy as np
from ast import literal_eval
from pykep.trajopt import mga_1dsm

import _JNSQ_System
import _Vanilla_System as _Vanilla_System

from _Kraken_Patch import transx 
from _Kraken_Patch import decode_dV_tof
from _Kraken_Patch import EJ_BurnDV
from _Kraken_Patch import EJ_Pe_Direction 
from _Kraken_Patch import EJ_angle_from_Pe
from _Kraken_Patch import EJ_angle_to_Prograde 
from _Kraken_Patch import EJ_details


from math import sqrt, pi, cos, sin, acos,log
from itertools import product
from collections import OrderedDict
import re, ast
from matplotlib import pyplot as plt, cm
import matplotlib.colors as colors
import seaborn as sns




'''
class mga_gene:
    #wrapper of the naked mga gene to avoid upsetting pandas
    def __init__(self,gene):
        #self._gene = np.zeros(len(gene))
        self._gene = gene
'''



class Wayfinder:
    '''
    This class will allow to search/scan one or several gravity assists sequences
    for one or several T0 and/or ToF bins.
    
    The data is stored in .xlsx format. All the object except dictionaries are rebuilt from the data.
    
    The main functions are : add_mga_scans => allow to add a bunch of fly_by_sequences to scan.
    The fly_by_sequences are defined with a combinatorial method, including wildcard to allow skipping a step,
    e.g. ; [["Kerbin"],["Eve"],["Kerbin","Eve,"*"],["Moho"]] means from Kerbin to Moho, passing by Eve as a first
    step, then Either Kerbin, Eve or ignore the second step.
    
    the save/load method will allow to save/load the df. Usually those are not called explicitly, since edit and optimize will save by default,
    and __init__ will call load.
    
    Note : optimize is pretty dumb at this point and can't target a batch to optimize yet. Il will go through everything in a "TODO" state.
    
    Be aware that using the 'high' optimizsation level can get very time consuming, especially in large batches.
    
    
    TODO / IDEA list :
        Add a grid plot of sequence vs t0 with DV as color ?
        quipu style plot of transfers with color coded DV and planets passed as pearls, length between knots for the time ?
        /!\ Ejection altitude as a job parameter and patched in the fitness eval via the decorated problem.
        => Full EJ DV calc is too costly, for now simplified EJ_DV calc (underestimate the ejection inclination+ward on i at SOI to ~25°)
        => A trick would be to pre-calc a map from iSOI+Vinf => iEJ (estimate)
        => Will do for now, but will need to revisit.
        => A neat trick would be to toggle the full EJDV calc only when already deep in the optimisation. but how ? can't see a way.
        Add tests => done
        Maybe add a smart toggle in the decorated problem to trigger the use of the full ej calc when the solution gets "good enough"
        
        
    '''

    


    def __init__(self,datastore_name = 'WayFinder_KSP',planet_pack = "Vanilla"):
        
        self.planet_pack = planet_pack
        if self.planet_pack == "Vanilla":
            self._Edy2Kdy = 4
            self._datastore_name = "WayFinder_Vanilla"
            self._Body_abrev_dic = {
                "Moho"       : 'Mo',
                "Eve"        : 'E' ,
                "Kerbin"     : 'K' ,
                "Duna"       : 'D' ,
                "Dres"       : 'Dr',
                "Jool"       : 'J',
                "Eeloo"      : 'El',
                }
            self._fullname_dic = {
                'Moho'       : _Vanilla_System.Moho   ,
                'Eve'        : _Vanilla_System.Eve    ,
                'Kerbin'     : _Vanilla_System.Kerbin ,
                'Duna'       : _Vanilla_System.Duna   ,
                'Dres'       : _Vanilla_System.Dres   ,
                'Jool'       : _Vanilla_System.Jool   ,
                'Eeloo'      : _Vanilla_System.Eeloo   ,
                }
        elif self.planet_pack == "JNSQ":
            self._Edy2Kdy = 2
            self._datastore_name = "WayFinder_JNSQ"
            self._Body_abrev_dic = {
                "Moho"       : 'Mo',
                "Eve"        : 'E' ,
                "Kerbin"     : 'K' ,
                "Duna"       : 'D' ,
                "Edna"       : 'Ed',
                "Dres"       : 'Dr',
                "Jool"       : 'J',
                }    
            self._fullname_dic = {
                'Moho'       : _JNSQ_System.Moho   ,
                'Eve'        : _JNSQ_System.Eve    ,
                'Kerbin'     : _JNSQ_System.Kerbin ,
                'Duna'       : _JNSQ_System.Duna   ,
                'Edna'       : _JNSQ_System.Edna   ,
                'Dres'       : _JNSQ_System.Dres   ,
                'Jool'       : _JNSQ_System.Jool   ,
                }

        self._opt_levels_dic = {
                'debug'      : [1,1,1,1],             # 0   UT debuging mode
                'low'        : [100,7,25,140],         # 1   UT
                'moderate'   : [100,16,25,140],       # 3   UT
                'high'       : [120,24,25,160],       # 4.5 UT   
                'wide'       : [100,48,50,140],       # 6   UT
                'ultra'      : [150,48,25,200],       # 11  UT
                'ultra+'     : [150,64,35,200],       # 19  UT         
                'deep'       : [150,24,25,200],       # 5.5 UT
                }
        
        self._opt_insertion_dic = {
                # "add_vinf_arr", "mga_orbit_insertion", "mga_e_target", "mga_alt_target" 
                'circular'          : [True, True ,0.0],   #
                'elliptical'        : [True, True ,0.9],   #
                'vinf'              : [True ,False,0.0],   #
                'none'              : [False,False,0.0],   #   
                }
        
        
        '''Monkey patching mga_1dsm for vanilla/JNSQ'''  
        mga_1dsm.EJ_BurnDV              = EJ_BurnDV
        mga_1dsm.EJ_Pe_Direction        = EJ_Pe_Direction
        mga_1dsm.EJ_angle_from_Pe       = EJ_angle_from_Pe
        mga_1dsm.EJ_angle_to_Prograde   = EJ_angle_to_Prograde
        mga_1dsm.EJ_details             = EJ_details
        mga_1dsm.transx                 = transx    
        mga_1dsm.decode_dV_tof          = decode_dV_tof 

        
    def add_batch(self,swing_by_bodies,
                 t0_min=0  , t0_bin=100 ,  n_t0_bins=10,
                 tof_min=0 , tof_bin=1500, n_tof_bins=1,
                 opt_level = 'debug', datastore_name = 'default', 
                 auto_tof = False, 
                 overwrite = False,
                 opt_injection = 'circular',
                 injection_altitude = 100000, 
                 ejection_altitude = 100000,
                 lambert_max_revs = 0) :
        print(type(lambert_max_revs))
        
        
        ''' 
        First, check if there is an existing file
        if there is, then load it.
        '''
        if datastore_name != "default":
            self._datastore_name = datastore_name
        
        if os.path.isfile(self._datastore_name+'.xlsx'):
            print('There is an existing file with name : '+self._datastore_name+'.xlsx')
            self.load_df(self._datastore_name)

        ''' T0s and tofs don't need to be part of the object.'''
        T0s  = range(t0_min,t0_min+n_t0_bins*t0_bin,t0_bin)
        ToFs = range(tof_min,tof_min+n_tof_bins*tof_bin,tof_bin)     
        
        sequences = self.generateSequences(swing_by_bodies) 
        shortSequences = self.generateShortSequences(swing_by_bodies)
        
        seqs2seqnames_dic = dict(zip(shortSequences, sequences))
        if auto_tof:
            idx_tuples = []
            for seq_shortname in shortSequences:
                tof_min,tof_bin     = self.auto_tof(seqs2seqnames_dic[seq_shortname])
                n_tof_bins  = 1
                for t0 in T0s:
                    idx_tuples.append((seq_shortname,t0,tof_min))
            idx = pd.MultiIndex.from_tuples(idx_tuples, names=['Seq', 'T0_lb', 'ToF_lb'])

        else :
            idx = pd.MultiIndex.from_product([shortSequences,T0s,ToFs],
                                          names=['Seq', 'T0_lb', 'ToF_lb'])          
        
        
        
        col = ['gene',                                                          # Output gene
               'job_status','job_sade_gen','job_n_island','job_island_pop','job_n_evo_steps', # Job informations
               'mga_seq_shortname','mga_seq_fullname','mga_tof','mga_t0','mga_vinf','mga_tof_encoding',      # mga informations
               'mga_add_vinf_dep','mga_add_vinf_arr','mga_multi_objective',     #
               'mga_alt_start','mga_alt_target',                                # added ejection orbit altitude (1.2)
               'mga_orbit_insertion','mga_rp_target','mga_e_target',            #
               'mga_problem',                                                   # udp problem of type MGA (object)
               'batch_t0_min','batch_t0_bin','batch_n_t0_bins',                 # batch t0  information
               'batch_tof_min','batch_tof_bin','batch_n_tof_bins',              # batch tof information
               'batch_sequences','batch_opt_level']                             # batch sequence
        
        '''All the mga problem parameters are stored flat in the DF. sequence is stored as a 
        list of string and require conversion when instancing'''

        add_df = pd.DataFrame('-',idx, col)

        
        add_df.astype('object')
        #now that the DF is set, put the optimization parameters in
        add_df.loc[:,'job_sade_gen']    = self._opt_levels_dic[opt_level][0]
        add_df.loc[:,'job_n_island']    = self._opt_levels_dic[opt_level][1]
        add_df.loc[:,'job_island_pop']  = self._opt_levels_dic[opt_level][2]
        add_df.loc[:,'job_n_evo_steps'] = self._opt_levels_dic[opt_level][3]
        add_df.loc[:,'job_status']      = 'TODO'


        '''
        in order to store things in a 'safe'& clear fashion, the best is to limit onself to data and avoid objects.
        This means that we need to store genes and sequences as lists.
        
        Filling the df using a loop is likely a bad idea performance-wise but... we'll see.
        '''
        print(add_df)
        print("****************************************************")
        
        for seq_shortname in shortSequences :
            if auto_tof:
                tof_min,tof_bin     = self.auto_tof(seqs2seqnames_dic[seq_shortname])
                n_tof_bins  = 1  
                
            ToFs        = range(tof_min,tof_min+n_tof_bins*tof_bin,tof_bin)
            target      = self._fullname_dic[seqs2seqnames_dic[seq_shortname][-1]]  
                                
            for t0 in T0s :
                for tof in ToFs:
                    add_df.at[(seq_shortname,t0,tof),'batch_t0_min']        = t0_min                   
                    add_df.at[(seq_shortname,t0,tof),'batch_t0_bin']        = t0_bin  
                    add_df.at[(seq_shortname,t0,tof),'batch_n_t0_bins']     = n_t0_bins  
                    add_df.at[(seq_shortname,t0,tof),'batch_tof_min']       = tof_min  
                    add_df.at[(seq_shortname,t0,tof),'batch_tof_bin']       = tof_bin  
                    add_df.at[(seq_shortname,t0,tof),'batch_n_tof_bins']    = n_tof_bins    
                    add_df.at[(seq_shortname,t0,tof),'batch_sequences']     = swing_by_bodies
                    add_df.at[(seq_shortname,t0,tof),'batch_opt_level']     = opt_level
                    add_df.at[(seq_shortname,t0,tof),'batch_opt_insertion'] = opt_injection
                    add_df.at[(seq_shortname,t0,tof),'mga_seq_shortname']   = seq_shortname
                    add_df.at[(seq_shortname,t0,tof),'mga_seq_fullname']    = seqs2seqnames_dic[seq_shortname]
                    add_df.at[(seq_shortname,t0,tof),'mga_t0']              = copy.deepcopy([t0/self._Edy2Kdy,(t0+t0_bin)/self._Edy2Kdy])
                    add_df.at[(seq_shortname,t0,tof),'mga_tof']             = copy.deepcopy([tof/self._Edy2Kdy,(tof+tof_bin)/self._Edy2Kdy])
                    add_df.at[(seq_shortname,t0,tof),'mga_vinf']            = [0.8, 1.8]
                    add_df.at[(seq_shortname,t0,tof),'mga_tof_encoding']    = 'alpha'
                    add_df.at[(seq_shortname,t0,tof),'mga_add_vinf_dep']    = True                                  #This should always be true, can't see a use case without it
                    add_df.at[(seq_shortname,t0,tof),'mga_add_vinf_arr']    = self._opt_insertion_dic[opt_injection][0]
                    add_df.at[(seq_shortname,t0,tof),'mga_multi_objective'] = False                                 #Setting to True will mess things up, not useable for now.
                    add_df.at[(seq_shortname,t0,tof),'mga_orbit_insertion'] = self._opt_insertion_dic[opt_injection][1]                                  
                    add_df.at[(seq_shortname,t0,tof),'mga_alt_start']       = ejection_altitude                     #NB : not part of mga definition per se, but gets into the decorator
                    add_df.at[(seq_shortname,t0,tof),'mga_alt_target']      = injection_altitude                    #NB : only for bookeeping and help readability
                    add_df.at[(seq_shortname,t0,tof),'mga_rp_target']       = self.rp_target_ward(target,injection_altitude)    #NB : rp target is radius+parking alt of target. 
                    add_df.at[(seq_shortname,t0,tof),'mga_e_target']        = self._opt_insertion_dic[opt_injection][2]
                    add_df.at[(seq_shortname,t0,tof),'mga_lambert_max_revs'] = int(lambert_max_revs)
                    add_df.at[(seq_shortname,t0,tof),'gene']                = []  
                    add_df.at[(seq_shortname,t0,tof),'result_DV']           = 99999
                    add_df.at[(seq_shortname,t0,tof),'result_t0']           = 99999
                    add_df.at[(seq_shortname,t0,tof),'result_tof']          = 99999               
                    add_df.at[(seq_shortname,t0,tof),'result_ej_vinf']      = 99999 

                    #self.build_one_mgaproblem(seq_name,t0,tof) #this shoud be done after the merge part.
                    
        
        '''
        Now that we have the add_df frame, we can either assign it or merge it to the self._df
        '''
        print(add_df)
        print("****************************************************")
        if os.path.isfile(self._datastore_name+'.xlsx'):
            self.load_df(self._datastore_name)
            frames = [self._df, add_df]            
            self._df = pd.concat(frames)
            
            '''
            We remove duplicates in terms of row index
            If the 'overwrite' option is on, the added line will erase existing ones, 
            including optimisation results...
            '''
            
            if overwrite :
                self._df = self._df[~self._df.index.duplicated(keep='last')]                
            #else :
            #    self._df = self._df[~self._df.index.duplicated(keep='first')]
            
        else : 
            self._df = add_df
            
        print(self._df)
        
        self.save_df() #save the changes.
            
         
            
        
        #print(self._df.dtypes)
        
        #self.build_all_mgaproblems()


    def edit_batch(self,swing_by_bodies,action = 'reset',save_it = True):
            '''
            v1.1.1
            New function to edit all jobs for a given "swingby" list. 
            should works as follow: takes a swingby sequence and an option :
                - 'reset' => changes the job status to TODO
                - 'debug'/'low'/'moderate'/'high' => change the opt level and reset to TODO
            
            Later / nice to have : ability to split/merge bin of a batch by a factor of two.
            This is more complicated as it means : droping the df entries, and adding the new ones.
            
            low priority : ability to use shorthand notation when adressing a single sequence (e.g. KED instead of [['Kerbin']['Eve']['Duna']])
            => do in the same way it is done for the other methods.
            '''
            
            ''' We generate the sequences from the swing_by_bodies'''
            sequences = self.generateSequences(swing_by_bodies) 
            '''note : This is messy and there must be a smarter way => done with separator and split'''
            for seq in sequences :
                for MI, new_df in self._df.groupby(level=[0,1,2]):
                    #print(MI)
                    if self._df.loc[MI,'mga_seq_fullname'] == seq:
                        if action == 'reset':
                            self._df.at[MI,'job_status']        = 'TODO'                        
                        elif 'lvl:' in action :
                            self._df.at[MI,'batch_opt_level'] = action.split(':')[-1]
                            print("Optimization level set to " +action.split(':')[-1])
                            self._df.at[MI,'job_status']        = 'TODO'                                                                
                        elif 'inj:' in action :
                            self._df.at[MI,'batch_opt_insertion'] = action.split(':')[-1]
                            print("Injection type set to " +action.split(':')[-1])
                            self._df.at[MI,'job_status']        = 'TODO'    
                        else : 
                            print("Error : did not recognised action")                        
                        
                        
                        self._df.at[MI,'job_sade_gen']          = self._opt_levels_dic[self._df.at[MI,'batch_opt_level']][0]
                        self._df.at[MI,'job_n_island']          = self._opt_levels_dic[self._df.at[MI,'batch_opt_level']][1]
                        self._df.at[MI,'job_island_pop']        = self._opt_levels_dic[self._df.at[MI,'batch_opt_level']][2]
                        self._df.at[MI,'job_n_evo_steps']       = self._opt_levels_dic[self._df.at[MI,'batch_opt_level']][3]
                        self._df.at[MI,'mga_add_vinf_arr']      = self._opt_insertion_dic[self._df.at[MI,'batch_opt_insertion']][0]
                        self._df.at[MI,'mga_orbit_insertion']   = self._opt_insertion_dic[self._df.at[MI,'batch_opt_insertion']][1]
                        self._df.at[MI,'mga_e_target']          = self._opt_insertion_dic[self._df.at[MI,'batch_opt_insertion']][2]                            
                        #print(self._df.loc[MI,'job_status'])
                        #print(self._df.loc[MI,'batch_opt_level'])
                        
            if save_it:
                self.save_df()


    def build_one_mgaproblem(self,seq_name,t0,tof,):
        ''' 
        This method re-creates the mga problem from the data using a dict mapping the planet names to the planets objects.    
        Does it for one item (seq,t0,tof)
        '''
        
        
        
        self._df.loc[(seq_name,t0,tof),'mga_problem']         = copy.deepcopy(mga_1dsm(
                            seq                     = list(map(self._fullname_dic.get, self._df.loc[(seq_name,t0,tof),'mga_seq_fullname'])),
                            t0                      = copy.deepcopy(self._df.loc[(seq_name,t0,tof),'mga_t0']),
                            tof                     = copy.deepcopy(self._df.loc[(seq_name,t0,tof),'mga_tof']),
                            vinf                    = self._df.loc[(seq_name,t0,tof),'mga_vinf'],
                            tof_encoding            = self._df.loc[(seq_name,t0,tof),'mga_tof_encoding'],
                            add_vinf_dep            = self._df.loc[(seq_name,t0,tof),'mga_add_vinf_dep'],
                            add_vinf_arr            = self._df.loc[(seq_name,t0,tof),'mga_add_vinf_arr'],
                            multi_objective         = bool(self._df.loc[(seq_name,t0,tof),'mga_multi_objective']),
                            orbit_insertion         = self._df.loc[(seq_name,t0,tof),'mga_orbit_insertion'],
                            rp_target               = self._df.loc[(seq_name,t0,tof),'mga_rp_target'],
                            e_target                = self._df.loc[(seq_name,t0,tof),'mga_e_target'],
                            max_revs                = int(self._df.loc[(seq_name,t0,tof),'mga_lambert_max_revs'])))
        
    
    def build_all_mgaproblems(self):
        ''' 
        This method re-creates the mga problem from the data using a dict mapping the planet names to the planets objects.      
        Does it for all rows in the DF
        '''

        for MI, new_df in self._df.groupby(level=[0,1,2]):
            self._df.at[MI,'mga_problem']          = copy.deepcopy(mga_1dsm(
                            seq                     = list(map(self._fullname_dic.get, self._df.loc[MI,'mga_seq_fullname'])),
                            t0                      = copy.deepcopy(self._df.loc[MI,'mga_t0']),
                            tof                     = copy.deepcopy(self._df.loc[MI,'mga_tof']),
                            vinf                    = self._df.loc[MI,'mga_vinf'],
                            tof_encoding            = self._df.loc[MI,'mga_tof_encoding'],
                            add_vinf_dep            = self._df.loc[MI,'mga_add_vinf_dep'],
                            add_vinf_arr            = self._df.loc[MI,'mga_add_vinf_arr'],
                            multi_objective         = bool(self._df.loc[MI,'mga_multi_objective']), 
                            orbit_insertion         = self._df.loc[MI,'mga_orbit_insertion'],
                            rp_target               = self._df.loc[MI,'mga_rp_target'],
                            e_target                = self._df.loc[MI,'mga_e_target'],
                            max_revs                = int(self._df.loc[MI,'mga_lambert_max_revs'])))


        
    def save_df(self):
        print(self._datastore_name)
        self._df.to_excel(self._datastore_name+".xlsx")
        
    def load_df(self,datastore_name = "Default"):
        if datastore_name == "Default":
            if self.planet_pack == "Vanilla":
                datastore_name = "WayFinder_Vanilla"
            elif self.planet_pack == "JNSQ":
                datastore_name = "WayFinder_JNSQ"
        '''
        #Here we need to finish building the full MGA_scan object from the data in the df.
        #This means we need to create the proper mga object for each entry.
        '''                     
        
        
        
        self._datastore_name = datastore_name
        self._df = pd.read_excel(self._datastore_name+".xlsx",index_col=[0,1,2]) 
        self._df['mga_seq_fullname'] = copy.deepcopy(self._df['mga_seq_fullname'].apply(literal_eval))
        self._df['mga_t0'] = self._df['mga_t0'].apply(literal_eval)
        self._df['mga_tof'] = self._df['mga_tof'].apply(literal_eval)
        self._df['mga_vinf'] = self._df['mga_vinf'].apply(literal_eval)
        self._df['gene'] = self._df['gene'].apply(literal_eval)
        #self._df['gene'] = pd.eval(self._df['gene'])    
        
        self.build_all_mgaproblems()
        


    def debugPrint(self):
        for MI, new_df in self._df.groupby(level=[0,1,2]):
            print(self._df.loc[MI,'mga_seq_fullname'])

    def rp_target_ward(self,target,injection_altitude):
        '''This function is a ward to avoid injection below the safe altitude or in atmosphere'''
        '''It will display a warning, and set a higher injection altitude of 1000 km'''
        if target.radius+injection_altitude < target.safe_radius:            
            print("Warning : injection altitude below safe limits, injection altitude for target ("+target.name+") set to 1000 km")
            return target.radius+1000000
        else :
            return target.radius+injection_altitude

    def generateSequences(self,swing_by_bodies=[["Kerbin"],["Eve","Duna"],["Kerbin"]]):
        '''
        Generates flyby sequences possibilities from start to end passing by the different combinations of bodies in the swing_by_bodies list
        The overall sequence length will be at minimum 3 and at most equal to the length parameter
        '''
        
        seqs = list(i for i in product(*swing_by_bodies) if ((len(set(i[1:-1])) >= len(i[1:-1])-1 or len(i)-i.count("*")) < 5 and len(set(i[0:3])) > 0))
        
            
        sequences = []
        for i in seqs:
            temp = []
            for j in i:
                if isinstance(j, list):
                    for k in j:
                        temp.append(k)
                else:
                    temp.append(j)
            sequences.append(temp)

        '''Cleanup the *'s              '''
        
        for sbb in range(len(swing_by_bodies)-2):
            for i in range(len(sequences)):
                try:
                    sequences[i].remove('*')
                except ValueError:
                    pass  # do nothing!

        '''Remove duplicates            '''
        
        new_k = []
        for elem in sequences:
            if elem not in new_k:
                new_k.append(elem)
        sequences = new_k
        
        return sequences
        
    def generateShortSequences(self,swing_by_bodies=[["Kerbin"],["Eve","Duna"],["Kerbin"]]):
        '''
        Generate the list in shorthand notation for the flyby sequences.      
        '''
        shortSequences = []
        sequences = self.generateSequences(swing_by_bodies)
        for seq in sequences :
            seq_name = ''
            for body in seq:
                #print(body)
                seq_name += self._Body_abrev_dic[body]
            shortSequences.append(seq_name)
        return shortSequences
    
    def orbital_period(self,planet):
        JNSQ_Dy2s = 12*3600
        Vanilla_Dy2s = 6*3600
        period = 2 * pi * sqrt( planet.orbital_elements[0]**3 / planet.mu_central_body )
        if self.planet_pack == "JNSQ":
            return round(period / JNSQ_Dy2s,0)
        elif self.planet_pack == "Vanilla":
            return round(period / Vanilla_Dy2s,0)            
        else :
            print("Unknown planet pack, auto-tof will use Vanilla ksp values for its guess !")
            return round(period / Vanilla_Dy2s,0)  
        
    def auto_tof(self,seq_fullname,debug = False):
        '''
        Improving the auto_tof feature
        1) when the sequence contains KXXT, with XX a dual assist, plan 1-2Xyr margin.
        2) For the last body, it should never be a full year, half at most, so 0.2-0.4 of period
        '''
        planet_sequence = list(map(self._fullname_dic.get, seq_fullname))
        tof_guess = 0

        for prev_planet, planet in zip(planet_sequence[0:], planet_sequence[1:]):

            if planet == prev_planet:
                tof_guess += 2*self.orbital_period(planet)
            elif planet == planet_sequence[-1] :
                tof_guess += 0.4*self.orbital_period(planet)
            elif planet != prev_planet:  
                tof_guess += self.orbital_period(planet)
            if debug:
                print(str(planet.name)+" "+str(tof_guess))

        '''
        for planet in planet_sequence[1:]:
            tof_guess += self.orbital_period(planet)
            if debug:
                print(str(planet.name)+" "+str(self.orbital_period(planet)))
        '''
        '''I chose tof between 0.55 and 1.1 of the sum'''
        if debug :
            print("tof_lb : "+str(int(round(tof_guess*0.5,-2)))+" tof_ub : "+str(int(round(tof_guess*1.0,-2))))
        return int(round(tof_guess*0.5,-2)), int(round(tof_guess*1.0,-2))
   
    def recalc_results(self):
                
        '''
        This methid will just recalculate the DV, T0 and tof and ejection Vinf of the result and store the result.
        '''
        for MI, new_df in self._df.groupby(level=[0,1,2]):
            if self._df.loc[MI,'job_status'] == 'DONE':
                udp = self._df.at[MI,'mga_problem']
                self._df.at[MI,'result_DV'],self._df.at[MI,'result_t0'],self._df.at[MI,'result_tof'],self._df.at[MI,'result_ej_vinf'] = udp.decode_dV_tof(self._df.at[MI,'gene'])     
        self.save_df()
        
    def audit_results(self):
        '''
        This method will look at the TOFs and Vinf and check if some completed jobs are hitting a boundary.
        '''
        for MI, new_df in self._df.groupby(level=[0,1,2]):
            if self._df.loc[MI,'job_status'] == 'DONE':
                if self._df.at[MI,'result_tof'] >= self._df.at[MI,'mga_tof'][1]*4:
                    print("upper tof boundary hit for job :"+str(MI))
                if self._df.at[MI,'result_tof'] <= self._df.at[MI,'mga_tof'][0]*4:
                    print("lower tof boundary hit for job :"+str(MI))
                if self._df.at[MI,'result_ej_vinf'] >= self._df.at[MI,'mga_vinf'][1]*1000:
                    print("upper ej vinf boundary hit for job :"+str(MI))
                if self._df.at[MI,'result_ej_vinf'] <= self._df.at[MI,'mga_vinf'][0]*1000:
                    print("lower ej vinf boundary hit for job :"+str(MI))                    
     
    def optimize(self, n=5000, save_it = True):
        count = 0
        '''
        We need to pull out the sequences we plan to optimize.
        Each batch has the sequence in a string format, so that should work.
        A problem is when there will be multiple batches in a single file.
        How to pick up ?
        
        Maybe for now just brute parse all the jobs, indifferent of the batch and 
        exectute everything that is in a TODO state.
        
        By default it will save the results at the end of the optimisation process
        '''

        ''' Parse on the multiindex using grouping, then run the optimisation '''
        
        for MI, new_df in self._df.groupby(level=[0,1,2]):
                                    
            if count < n and self._df.loc[MI,'job_status'] == 'TODO':
                udp           = self._df.loc[MI,'mga_problem']  
                sade_gen      = int(self._df.loc[MI,'job_sade_gen'])
                n_island      = int(self._df.loc[MI,'job_n_island'])
                island_pop    = int(self._df.loc[MI,'job_island_pop'])
                n_evo_steps   = int(self._df.loc[MI,'job_n_evo_steps'])
                mga_alt_start = int(self._df.loc[MI,'mga_alt_start'])
                '''There we set the problem to solve as the UDP'''
                
                '''v1.2: In order to improve accuracy, a decorated problem is used to include the ejection burn cost
                including oberth. For simplicity and speed, we assume a planar burn (zero inclination) as an approximation'''
                
                '''v1.2: For some reason, if the decorator is not defined within the optimization loop,
                the optimizer does not iterate past the first iteration. There must be a better/cleaner
                way to do it, but I cannot see it for now'''
                
                '''v1.2: we need to access the fullnamedict in the decorator, 
                to get the right planetary parameters, there must be a smarter way to do it but can't see it'''
                
                '''v1.3 ejection burn altitude is now a proper parameter, instead of a default 100km
                Also, "tweaks" to the decorated problem fitness : 
                    - KK inside a sequence will force 2yr period 
                    - Ejections with inclinations above 15° will are pruned
                    - Considered (to be benchmarked) : KK time constraint to n*1/2 Kyr if at the start, eg in KKEMo
                    
                    '''
                
                
                encap_fullnamedic = self._fullname_dic
                
                def f_decor_vanilla(orig_fitness_function):
                    def new_fitness_function(self, dv):

                        '''
                        This is patchy way of dealing with the problem of not accounting the full ejection burn
                        we remove the initial v_inf from the fitness and replace it by the ejction burn value, assuming it's a planar burn.
                        '''                          
                        sequence_rgx = re.compile("\[.*?\]")
                        bracket = sequence_rgx.search(self.inner_problem.get_extra_info())
                        sequence =  ast.literal_eval(bracket.group())
                        FirstBody = encap_fullnamedic[sequence[0]]
                        
                        theta = 2 * pi * dv[1]
                        phi = acos(2 * dv[2] - 1) - pi / 2
                        Vinfx = dv[3] * cos(phi) * cos(theta)
                        Vinfy = dv[3] * cos(phi) * sin(theta)
                        Vinfz = dv[3] * sin(phi)
                        Vinfxy = sqrt(Vinfy**2+Vinfx**2)
                        

                        Rsoi = FirstBody.orbital_elements[0] * (FirstBody.mu_self/FirstBody.mu_central_body)**(0.4)
                        v0 = sqrt(FirstBody.mu_self / (FirstBody.radius+mga_alt_start) )
                        vxy_ob = sqrt(Vinfxy**2 + 2 * FirstBody.mu_self / (FirstBody.radius+mga_alt_start) -2 * FirstBody.mu_self / Rsoi) - v0 #oberth only for xy
                        v_ob = sqrt(vxy_ob**2+Vinfz**2)
                        
                        if v_ob < (vxy_ob+abs(Vinfz)):
                            fitness = orig_fitness_function(self, dv)+v_ob-dv[3]
                        else :
                            fitness = orig_fitness_function(self, dv)+vxy_ob+abs(Vinfz)-dv[3]    

                        """prune solutions with SOI inclinations above 15°"""
                        fitness += max(0,abs(Vinfz/Vinfxy)-0.25)*10000
                        
                        """if the sequence contains a K-K slingshot this constraint prunes trajectories with KK tofs other than 2yr"""                        
                        """if KK is at the start, then ask for tof of 0.5 Kyr to 1.5 Kyr"""
                        KYr = 426/4.
                        HalfKyrWindow = KYr*0.52
                        CropWindow = 5/4.
                        T = list([0] * (len(sequence)-1))
                        for i in range(len(T)):
                            T[i] = -log(dv[5 + 4 * i])
                            alpha_sum = sum(T)

                        retval_T = [dv[-1] * time / alpha_sum for time in T]
                
                        if sequence[0] == "Kerbin" and sequence[0] == sequence[1]:
                            tof = retval_T[0]
                            #fitness += (max(CropWindow,min(tof%HalfKyr,abs(tof%HalfKyr-HalfKyr)))-CropWindow)*10000
                            fitness += 10000*(max(HalfKyrWindow,abs(tof-KYr))-HalfKyrWindow)

                        for planet, next_planet, tof in zip(sequence[1:], sequence[2:], retval_T[1:]):
                            if planet == "Kerbin" and next_planet == planet:
                                fitness += 10000*(max(CropWindow,abs(tof-KYr*2))-CropWindow)
                
                        return fitness
                    return new_fitness_function
                
                def f_decor_jnsq(orig_fitness_function):
                    def new_fitness_function(self, dv):

                        '''
                        This is patchy way of dealing with the problem of not accounting the full ejection burn
                        we remove the initial v_inf from the fitness and replace it by the ejction burn value, assuming it's a planar burn.
                        '''                          
                        sequence_rgx = re.compile("\[.*?\]")
                        bracket = sequence_rgx.search(self.inner_problem.get_extra_info())
                        sequence =  ast.literal_eval(bracket.group())
                        FirstBody = encap_fullnamedic[sequence[0]]
                        
                        theta = 2 * pi * dv[1]
                        phi = acos(2 * dv[2] - 1) - pi / 2
                        Vinfx = dv[3] * cos(phi) * cos(theta)
                        Vinfy = dv[3] * cos(phi) * sin(theta)
                        Vinfz = dv[3] * sin(phi)
                        Vinfxy = sqrt(Vinfy**2+Vinfx**2)
                        

                        Rsoi = FirstBody.orbital_elements[0] * (FirstBody.mu_self/FirstBody.mu_central_body)**(0.4)
                        v0 = sqrt(FirstBody.mu_self / (FirstBody.radius+mga_alt_start) ) 
                        vxy_ob = sqrt(Vinfxy**2 + 2 * FirstBody.mu_self / (FirstBody.radius+mga_alt_start) -2 * FirstBody.mu_self / Rsoi) - v0 #oberth only for xy
                        v_ob = sqrt(vxy_ob**2+Vinfz**2)
                        
                        if v_ob < (vxy_ob+abs(Vinfz)):
                            fitness = orig_fitness_function(self, dv)+v_ob-dv[3]
                        else :
                            fitness = orig_fitness_function(self, dv)+vxy_ob+abs(Vinfz)-dv[3]    
                        
                        
                        """prune solutions with SOI inclinations above 15°"""
                        fitness += max(0,abs(Vinfz/Vinfxy)-0.25)*10000
                        
                        """if the sequence contains a K-K slingshot this constraint prunes trajectories with KK tofs other than 2yr"""
                        KK2yr = 2*365/2
                        KKWindow = 5/2.
                        T = list([0] * (len(sequence)-1))
                        for i in range(len(T)):
                            T[i] = -log(dv[5 + 4 * i])
                            alpha_sum = sum(T)
                            retval_T = [dv[-1] * time / alpha_sum for time in T]
                            for planet, next_planet, tof in zip(sequence, sequence[1:], retval_T):
                                if planet == "Kerbin" and next_planet == planet:
                                    fitness += 10000*(max(KKWindow,abs(tof-KK2yr))-KKWindow)
                                    
                        return fitness
                    return new_fitness_function


                if self.planet_pack == "JNSQ" :
                    dudp = pg.problem(pg.decorator_problem(udp,fitness_decorator=f_decor_jnsq))
                elif self.planet_pack == "Vanilla" :
                    dudp = pg.problem(pg.decorator_problem(udp,fitness_decorator=f_decor_vanilla))
                print(dudp)
                #pg.problem(udp)
                '''We solve it!!'''
                uda = pg.sade(gen=sade_gen)
                archi = pg.archipelago(algo=uda, prob=dudp, n=n_island, pop_size=island_pop)
                print(
                        "Running Sade Algo on "+str(n_island)+" islands")
                archi.evolve(n_evo_steps)
                archi.wait()
                print("----------------------------------------------------------")
                sols = archi.get_champions_f()
                idx = sols.index(min(sols))
                print("Done!! Solutions found are: ", archi.get_champions_f())
                self._df.at[MI,'gene'] = copy.copy(list(archi.get_champions_x()[idx]))
                self._df.at[MI,'job_status'] = copy.copy('DONE')
                """Here we get the champion t0,tof and DV and store them in the df"""
                self._df.at[MI,'result_DV'],self._df.at[MI,'result_t0'],self._df.at[MI,'result_tof'],self._df.at[MI,'result_ej_vinf'] = udp.decode_dV_tof(self._df.at[MI,'gene'])
                count += 1
                print("iteration n°"+str(count)+" Seq : "+str(MI[0])+" lb_t0: "+str(MI[1])+" lb_tof: "+str(MI[2]))
        if save_it :
            self.save_df()


                    
    def decode_solutions(self,swing_by_bodies = []):       
        filtered_df = self.sequence_filter(swing_by_bodies)

        for MI, new_df in filtered_df:
            if self._df.loc[MI,'job_status'] == 'DONE':
                print("-----------------------------------------------------------------------")
                print("decoding solution with index "+str(MI))
                udp = self._df.at[MI,'mga_problem']
                gene = self._df.at[MI,'gene']
                alt = self._df.at[MI,'mga_alt_start']
                udp.transx(gene,alt)      



           
                    
    def find_best_plan(self,swing_by_bodies,t0_range = [0,10000]):
        '''1.1 added time window argument to select over a restricted range of launch dates'''
        
        filtered_df, _ = self.sequence_filter(swing_by_bodies,t0_range,gby_levels=[0,1,2])
        #print(filtered_df.indices[0])
        best_dV = 999999

        for MI, new_df in filtered_df:
            if self._df.loc[MI,'job_status'] == 'DONE':
                udp = self._df.loc[MI,'mga_problem']
                gene = self._df.loc[MI,'gene']
                if udp.decode_dV_tof(gene)[0] < best_dV :
                    best_dV = udp.decode_dV_tof(gene)[0]
                    best_gene = copy.deepcopy(gene)
                    best_udp = udp       
                                            
        print("best dV is : "+str(round(best_dV,1)) + " m/s")
        best_udp.transx(best_gene)        
        
        
    def plot_by_sequences(self,swing_by_bodies,t0_range = [0,10000]):
        '''
        This will plot the best candidate per sequence
        v1.1.1 : Since the DF is meant to store loads of sequences,
             we need to filter what we plot. This means asking for a swing_by_list, 
             or a sequence, or start-end pair with a start = Kerbin by default.
             For a start we'll do the swing_by_bodies list approach, and the add the
             other options in v1.1.2 
             
        NB : (1.4) there is very likely a much leaner and pleasing way to do this. 
    
        '''
        filtered_df, shortSequences = self.sequence_filter(swing_by_bodies,t0_range,gby_levels=[0])

        best_dVs = np.zeros((len(filtered_df),3))        
        SQ_dict = dict(zip(shortSequences, range(len(filtered_df))))

        
        for SQ, new_df_1 in filtered_df: 
            best_dV = 999999        
            for MI, new_df_2 in new_df_1.groupby(level=[1,2]):
                if self._df.loc[(SQ,MI[0],MI[1]),'job_status'] == 'DONE':
                    udp = self._df.loc[(SQ,MI[0],MI[1]),'mga_problem']
                    gene = self._df.loc[(SQ,MI[0],MI[1]),'gene']
                    if udp.decode_dV_tof(gene)[0] < best_dV :
                        best_dVs[SQ_dict[SQ]] = udp.decode_dV_tof(gene)[0:3]
                        best_dV = udp.decode_dV_tof(gene)[0]
                        #print(str(udp.decode_dV_tof(gene))+" "+str(SQ)+" "+str(best_dV))
                        #print(str(SQ)+" "+str(MI)+" "+str(new_df_2.loc[(SQ,MI[0],MI[1]),'mga_seq_fullname']))

        
        
        Seq_vs_dVs = dict(zip(shortSequences, best_dVs))
        ''' Order by DV '''
        Seq_vs_dVs = OrderedDict(sorted(Seq_vs_dVs.items(), key=lambda t: t[1][0]))
        fig, ax = plt.subplots( figsize=(9,6))
        
        palette = colors.Normalize(vmin = min([i[1] for i in Seq_vs_dVs.values()]), vmax= max([i[1] for i in Seq_vs_dVs.values()]))
        plt.bar(range(len(Seq_vs_dVs)), np.stack(list(Seq_vs_dVs.values()),axis = 0)[:,0],color=plt.get_cmap('coolwarm')(palette(np.stack(list(Seq_vs_dVs.values()),axis = 0)[:,1])), 
                edgecolor='black', linewidth=1, align='center')
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=palette)
        sm._A = []
        ax.set_ylabel('Total DV')
        ''' color bar code '''
        cbar = fig.colorbar(sm, ax=ax)
        if self.planet_pack == "JNSQ" :
            cbar.set_label('ToF in JNSQ Days', rotation=90)
        elif self.planet_pack == "Vanilla" :
            cbar.set_label('ToF in KSP Days', rotation=90)        
        ''' named ticks '''
        plt.xticks(range(len(Seq_vs_dVs)), list(Seq_vs_dVs.keys()), rotation='vertical')
        fig.tight_layout()
        plt.show()
        
        
    def plot_DVvsT0(self,swing_by_bodies,t0_range = [0,10000]):
        print("displays DV vs T0 bins")
        """
        For all sequence, plot a line chart with DV vs T0, with graduations equal to the jobs bining
        ?? different bining => just pick the first and be done with it ? or 100 dy as default ?
        
        prefered way is to plot on a df directly : df.plot(color=df.columns, figsize=(5, 3))
        or with seaborn : sns.lineplot(data=df, x='x', y='y', hue='color')
        
        """
        
        filtered_df, shortSequences = self.sequence_filter(swing_by_bodies,t0_range,gby_levels=[0,1])        
        filtered_df = (self._df.loc[filtered_df["result_DV"].idxmin()])
        
        """The two commented lines below are ways to plot that did not work out but might be usefull an other time"""
        #filtered_df.pivot('result_t0', 'mga_seq_shortname', 'result_DV').interpolate(method='linear').plot()
        #filtered_df.groupby(level=[0]).plot(x="result_t0",y="result_DV")
        
        """In the end sns with (x,y,hue) did the best job and in a pleasing way"""
        plt.figure(figsize=(10, 8))
        sns.lineplot(x='result_t0', y='result_DV', hue='mga_seq_shortname', data=filtered_df)
        #plt.yscale('log')
        plt.legend(title='Sequence', bbox_to_anchor=(1.05, 1), loc='upper left')

        
        
    def sequence_filter(self,swing_by_bodies,t0_range = [0,10000],gby_levels=[0,1,2]):
        
        
        '''
        Three cases of input to deal with :
        1) a list of list of bodies sequence (combinatorial form)
        2) a list of sequence short hands
        3) a single sequence short hand
            
        t0_range is here to select a restricted range of launch dates.
        '''
        
        if isinstance(swing_by_bodies, list):
            if isinstance(swing_by_bodies[0], list):
                ''' Get the abreviated sequence from the swing_by_bodies list of lists, then filter'''
                shortSequences = self.generateShortSequences(swing_by_bodies)
                filtered_df = self._df[self._df["mga_seq_shortname"].isin(shortSequences)].query("T0_lb >= @t0_range[0] and T0_lb+batch_t0_bin <= @t0_range[1]").groupby(level=gby_levels)
            elif isinstance(swing_by_bodies[0], str):
                ''' This is already a short hand sequence list, so filter away...'''
                filtered_df = self._df[self._df["mga_seq_shortname"].isin(swing_by_bodies)].query("T0_lb >= @t0_range[0] and T0_lb+batch_t0_bin <= @t0_range[1]").groupby(level=gby_levels)
                shortSequences = swing_by_bodies
        elif isinstance(swing_by_bodies, str):
            ''' This is a single short hand sequence, so turn into list and filter away...'''
            swing_by_bodies = [swing_by_bodies]
            filtered_df = self._df[self._df["mga_seq_shortname"].isin(swing_by_bodies)].query("T0_lb >= @t0_range[0] and T0_lb+batch_t0_bin <= @t0_range[1]").groupby(level=gby_levels)
            shortSequences = [swing_by_bodies]            
        else :
            ''' If no sequence is provided, default to return all of it'''
            filtered_df = self._df.groupby(level=gby_levels)
            shortSequences = list(self._df.index.get_level_values(0).unique()) 

        return filtered_df, shortSequences