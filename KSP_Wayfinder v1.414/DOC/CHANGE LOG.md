
### v1.0 :
- working prototype
- data saved as pickle
    
        
### v1.1 :
    
- transition to xlsx for storage to 1) avoid safety issues 2) make the saved data "legible"
- Combinatorial method for sequence generation 
- reset methods for optimisation
- possibility to add new jobs to and existing file
- added edit method to allow re-reunning a job with a different opt level
- various code cleanup

### v1.2 :

- Added plot function to compare sequences
- default set of injection parameters (circular, highly elliptical, Vinfinite or none)
- injection and ejection altitudes are now fully fledged job parameters, and are accepted in the edit method too.
- auto Tof mode added
- added * in sequence for optional swing by
- added proper Ejection DV calculation to display on the flight plan.
	   
        
### v1.3 :

- optimization now works with a decorated MGA problem to help tweak the optimisation, for example : added an ejection burn estimation in the fitness function,
    will force a 2Kyr period for KK sequences except if the sequence is starting with KK. Experimental MoMo braking sequence helper. Curation of high inclination ejections.
- replaced the two kraken patches by a single one with planet pack toggle.
- added a T0 window selector in the sequence filter method
- added on \epsilon0 to ward off some rare errors when computing ejection angles and inclination (minimal precision impact)
- added test with a full ejection burn calculation in the optimisation function and proof it with direct transfers (success)
- improved transX style flight plan display, removal of useless parameters.
            
### v1.4 :

- when a job completes, t0, tof and DV values are stored in the file. This will make plotting much easier
- added DV vs T0 plot for aritrary sequences in all possible forms, including combinatorial. Very nice.
- added a few tests

### v1.414

- added examples of use for different cases
- re-structures folders to make things neatier
- smarter auto_tof. Will use 0.5-1 orbital per body in the sequence, except for :
	1) the starting one (not counted)
	2) if "planet-planet" is found in the sequence (but not the first two), 1.0-2.0 "planet" periods are added
	3) for the last body, only 0.2-0.4 period is counted.


