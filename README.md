# Wayfinder
Wayfinder is a Multiple Gravity Assist search tool for KSP. 

It's goal is to allow efficient search of gravity assist sequences using state of the art tools using pykep/pygmo
packages from ESA. At the moment it's a set of python scripts with no GUI.

## References and requirements : 
This work would not have been possible without 
- pykep : https://esa.github.io/pykep/ 
- pygmo : https://esa.github.io/pygmo2/
- transfer planer for JNSQ : https://github.com/LouisB3/ksp-lwp-jnsq
- the original transfer planer : https://github.com/alexmoon/ksp

Special thanks to ESA to make those package available to the public, it's awesome.

## How does it work :

Wayfinder uses a job batch system with the results of said jobs saved in an xslx format for storage and readability. 
Jobs can be added with the desired parameters (fly-by sequence, insertion type, search space bining, optimizsation level and so on).
Once added, jobs can be run in batches, and results will be saved. Once results are saved, the results can be searched and accessed with
a few utilites as : 

- finding the best result in a given set (sequence + launch dates)
- display several sequences for comparison in terms of DV cost and Time of flight.
- compare different sequences in a DV cost vs launch date line plot
- display an orrery style plot with the trajectory arcs
- display a job result as a flight plan in a text format

## Structure of the folders

- **Examples** provides a set 4 examples to get started and get a feel of what is possible and how. I recommend having a look there first :)
- **Wayfinder Vanilla** contains a library of optimized MGA's for Vanilla KSP, along with a benchmarking script and the main script used to build / use 
	the MGA library.
- **Wayfinder JNSQ** contains a (small) library of optimized MGA's for the JNSQ planetpack, along with a benchmarking script and the main script used to build / use 
	the MGA library.
- **Wayfinder core** contains the core of the code


## Setting up a job :

The amount of work required to setup decent jobs has been reduced to a minimum and streamlined as much as possible.
This does not mean that one cannot go beyond what is preset, but ease of use was the focus. As an example, here is the code to setup
a search for transfers to Moho using either a KEEMo or KEMo flyby sequence for lauch dates from 0 to 1000 with a search binning of 100 days :

	plans.add_batch(
            swing_by_bodies = [["Kerbin"],["Eve"],["Eve","*"],["Moho"]],
            t0_min          = 0,
            t0_bin          = 100,
            n_t0_bins       = 10,          
            auto_tof        = True,
            opt_level       = "moderate",
            opt_injection   = "circular",
           	injection_altitude = 100000)


Here, several options require explanations :

- **auto_tof** 			: option that will ask the program to compute a reasonnable guess for the time of flight boundaries, with a single bin.
- **opt_level** 		: option giving preset values to the optimizer for the search. 
	+ "debug" will not do anything usefull except to check that the job will run.
	+ "low" or "moderate" level will often be enough to explore possible solution. A moderate optimization can deal with simple sequences (3-4 bodies)
	+ "high" optimization level usually is able to deal apply with up to 5 bodies in a sequence.
	+ "wide" can help with exploration when the sequence becomes complex. 
	+ "ultra" will tax the CPU a lot, but can help when searching with wide windows and complex sequence. It's usually smarter to get a feel of the search space with simpler jobs first before moving to it.
- **opt_injection** 	: tell the program what type of injection it should work with during the optimization. options are : 
	+	"none", which means that the optimization will not care about injection cost in DV. 
	+	"vinf" means that only the velocity at infinity (SOI) will be counted. 
	+	"elliptical" will set parameters for injection into a highly elliptical orbit with eccentricity of 0.9. 
	+	"circular" injections means that the optimizer will include the cost of an injection into a circular orbit.
- **injection alt** 	: target altitude for the orbital injection at destination, if the orbit is circular or elliptical (pe).
- **ejection alt**  	: parking orbit altitude from which ejection will be done.
- **overwrite**			: if "True" configured jobs that have the same keys (sequence, min_t0 and min_tof) are overwritten. Otherwise not. Default is set to "False"
- **datastore_name** 	: if working with non default datastore (xlsx sheet), provide the name of the file.


Each job corresponds to a gravity assist sequence with boundaries on ToF and launch dates. Wayfinder can generate such jobs in batches,
and is able to guess reasonnable ToF boundaries if wished.

## Running a job :

Once a job has been configured, starting the optimizer is just a matter of loading the file and calling the 

	swing_by_bodies = [["Kerbin"],["Duna"],["Eve"],["Kerbin"]]
	plans.load_df()  
	plans.optimize(n=20)
	plans.find_best_plan(swing_by_bodies,t0_range=[0,4000])


Once the jobs are setup, Wayfinder can search for solutions for each job. Several optimization levels with varying computational costs are pre-configured. 
Job can be edited in batch to change the optimization level or the orbital insertion  if required.

Once optimized, the results (gene encoding the solution) are saved along with all the job parameters and the resulting total DV cost, time of flight and lauch date.

Wayfinder then allows to search the job library for the best route in a given lauch window and/or compare different possible routes (sequences) 

Trajectories are multiple gravity assists with deep space manoeuvres. There might be some differences and inaccuracies between KSP and Wayfinder, 
but this is usually a matter of tweaking manoeuvre nodes. While it can find very efficient routes, it may also make use of very small windows 
that may be tricky to nail down in KSP. In particular, some flyby's can occur with inclined orbits that can make finding the right encounter difficult. 
This is especially the case if the encounter occurs far in the future.

## Things that are hard :
	
Some sequences are harder than others. To alleviate (as in reduce the optimization time) some additional constraints are applied to the pygmo
problem by decorating it and patching the fitness. At the moment the following constraints have been added :
	
- Replacement of ejection Vinf by an approximation of the true ejection cost. To make it short, Vinf is not what we want to optimize on, as we seek
	the lowest trajectory from the parking orbit, and not from the starting body SOI. We assume a zero inclination parking orbit, and as such, only the planar components
	of the ejection Vinf will benefit from the Oberth effect. The normal componant on the other hand will increase when calculated in the parking orbit since the ejection
	vector from the parking orbit is not collinear to the ejection Vinf. A full calculation being too costly, an approximation is used.
- Kerbin-Kerbin slings are restrited to periods of 0.5 to 1.5 Kyr if the KK sling is at the start of the sequence (KKX) and two Kyrs if the KK slingshot is within the sequence (for example KEKKJ).

## Using a result to plan in KSP

The last part is using the found solution and replicating in KSP. While it seems simple, it's not always the case. Here are a few tips
and tricks :

- **Tight Windows:** some solutions have very thigth windows for encounters. This is especially the case for segments where the spaceship trajectory is not coplanar with the next body in the gravity assist chain. It might seems at time that there is no solution, or that the encounter is impossible, but it is often just a case of the enounter being (really) hard to pull off. Moho is the worst in that regard considering the high inclination differentials and tiny SOI.
- **Error Propagation:** when dealing with long chains of assists, it is often pointless to nail down the encounters past the next two. The reason is that each burn will require to readjust anyway due to the very high sensitivity to even 0.1m/s difference in a burn.
- **Deep Space Manoeuvres (DSM):** if the first window is not hit perfectly, the DSM values will start to differ from the actual required values. In short, it is important to try to get the encounters as right as possible and not be affraid to tweak the DSM's a bit.
- **Kerbin-Kerbin assists:** one cannot select kerbin as a target while in kerbin SOI. This makes those assist not trivial to pull off at times. The best way is to set a node at the time of the next encounter and tweak the trajectory using it as a guide.
- **Nodes as Time Stamps:** a very helpfull trick is to setup manoeuvre nodes as time stamps for the future encounters, as it help finding the setup for the different encounters and also help gauge how far from the prediction one is wandering. 
- **Beta Plane Angle** is given for each encounter. This is the angle between the plane containing the ship during an encounter and the ecliptic plane. This is massively helpfull to orient the flyby path correctly.
- **Pe/Ap of Arcs** are given. This can help a lot to see if the trajectory obtained is correct or not. If the Pe or Ap of a trajectory segment is wrong, then correction is required.




## Installation

I can only speak about my setup, which is winpython 3.6 and/or 3.7, along with pykep and pygmo.
