# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:56:18 2019

@author: v.fave
"""

"""
Defines JNSQ system in terms of pykep planets.

"""

#!/usr/bin/env ipython
# -*- coding: cp1252 -*-

import pykep as pk
import numpy as np
from collections import OrderedDict, namedtuple
from math import degrees, radians, pi, sqrt
from numpy import deg2rad
#from pykep import planet, epoch



"""
Kerbol system

In addition to representing the orbital elements through a numpy array,
representations are also created based on named tuples, and on PyKEP planets.
See: http://docs.python.org/library/collections.html#collections.namedtuple
See: http://keptoolbox.sourceforge.net/documentation.html#PyKEP.planet


USAGE EXAMPLES ('body' matrix):
	
	# accessing a column of 'body' by name
	body[:,col['a (km)']]	
	
	# accessing a body's orbital elements, by name
	body[ row['callisto'] ]
	
	# accessing a body's specific orbital element, by name
	body[ row['callisto'], col['a (km)'] ]

USAGE EXAMPLES ('body_tuple' dictionary of named tuples):
	
	>>> body_tuple.keys()
	['jupiter', 'io', 'europa', 'ganymede', 'callisto']
	
	>>> body_tuple['io']
	orbital_elements(body='io', Epoch=58849.0, a=422029.68714001001, e=0.0043085246617730003, i=0.040115486869660003, Node=-79.640061742992003, w=37.991267683986997, M=286.85240405644998)
	
	>>> body_tuple['io'].Node
	-79.640061742992003

USAGE EXAMPLES ('body_obj' dictionary of PyKEP.planet objects):
	
	>>> degrees( body_obj['io'].orbital_elements[3] )
	-79.640061742992
	
	>>> body_obj['io'].eph( pk.epoch_from_string('2025-01-01 00:00:00.000') )
	((377861107.98154724, 186406866.4808699, 283715.7154025446),
	 (-7615.668664758434, 15593.895465832084, -3.281760584438344))	
	


"""
G = 6.674e-11 #08

# Values from KSP wiki
# 	Table 4: Keplerian orbit elements of the kerbol planets at t=0
# 	Table 5: Satellite physical constants
# 	Table 6: Other constants and conversions
body = np.array( [
	# UID   Epoch    a (km)             e           i (deg)       Node (deg)    w (deg)     M (rad)  R (km)      Atm_h     Rs (km)    mu (km^3/s^2)     mu_c           R SOI (km)  
	[  0.,  0,   14522400.000 ,        0.2,          7.0,           70.0,        15.0,      0.0,    650.0,        0.0,      660.0,     1.80034e13*G,  1.25719e20*G,       26571.730  ],
	[  1.,  0,   27131000.000 ,        0.01,         2.1,           15.0,        45.0,      5.7,   2050.0,       70.0,     2120.0,       8.645e14*G,  1.25719e20*G,       233566.973 ],    
	[  2.,  0,   37525647.8984324,     0.02,         0.0,            0.0,         0.0,      0.0,   1600.0,       85.0,     1690.0,     3.76157e14*G,  1.25719e20*G,       231587.700 ],   
	[  3.,  0,   57189100.000 ,        0.051,        0.06,         135.5,       345.0,      0.9,    800.0,       60.0,      870.0,     3.19734e13*G,  1.25719e20*G,       131664.417 ],    
	[  4.,  0,   94080000.000 ,        0.07,         3.0,           30.0,       310.0,      0.0,    260.0,        0.0,      270.0,     7.94632e11*G,  1.25719e20*G,        49408.693 ],
	[  5.,  0,   112687000.000 ,       0.145,        5.0,          280.0,        90.0,      3.9,    360.0,        0.0,      370.0,     2.28515e12*G,  1.25719e20*G,        90298.250 ],     
	[  6.,  0,   189765000.000 ,       0.05,         1.304,         52.0,        30.0,      0.6,  14000.0,      500.0,    14600.0,     2.99515e16*G,  1.25719e20*G,      6745640.285 ],    
	[  7.,  0,   359571000.000 ,       0.03,         1.7,           80.0,        75.0,      3.3,   8000.0,      500.0,     8600.0,     8.83969e15*G,  1.25719e20*G,      7845111.216 ],     
	[  8.,  0,   471171300.000 ,       0.26,         6.15,          50.0,       260.0,      3.54,   600.0,        0.0,      610.0,     7.93456e12*G,  1.25719e20*G,       621194.499 ],
	[  9.,  0,   527129000.000 ,       0.1,          4.0,          165.0,       175.0,      4.7,    450.0,        0.0,      460.0,     2.97546e12*G,  1.25719e20*G,       469437.943 ],     
	[ 10.,  0,   1712000000.000 ,      0.35,         20.0,          90.0,       150.0,      2.5,   3600.0,      100.0,     3800.0,     1.90430e15*G , 1.25719e20*G,     20213314.821 ],      
	], dtype=np.float64 )  

# Note : QA'd to match the JNSQ transfer planer by louis3b. All seems okay.
    
#Kerbin = pk.planet.keplerian(pk.epoch(0), (13599840256 ,0,0,0,0,3.14), 1.1723328e18, 3.5316000e12, 600000, 670000 , 'Kerbin')
#Duna   = pk.planet.keplerian(pk.epoch(0), (20726155264 ,0.051,deg2rad(0.06) ,0,deg2rad(135.5),3.14), 1.1723328e18, 3.0136321e11,320000, 370000 , 'Duna')
# Mapping from name of a body to number of row in which it is represented in the body matrix
row = OrderedDict( [
	( 'Moho'   , 0 ),
	( 'Eve'    , 1 ),
	( 'Kerbin' , 2 ),    
	( 'Duna'   , 3 ),  
	( 'Edna'   , 4 ),      
	( 'Dres'   , 5 ),      
	( 'Jool'   , 6 ),  
	( 'Lindor' , 7 ), 
   ( 'Eeloo'  , 8 ), 
   ( 'Hamek'  , 9 ),
   ( 'Nara'   , 10),   
	] )


    
# Mapping from name of a value to number of column in which it is represented in the body matrix
col = OrderedDict( [
	( 'Object UID'    , 0 ),
	( 'Epoch'         , 1 ),	# epoch
	( 'a (km)'        , 2 ),	# semi major axis
	( 'e'             , 3 ),	# eccentricity
	( 'i (deg)'       , 4 ),	# inclination
	( 'Node (deg)'    , 5 ),	# longitude of the ascending node
	( 'w (deg)'       , 6 ),	# argument of periapsis
	( 'M (rad)'       , 7 ),	# mean anomaly at epoch
	( 'R (km)'        , 8 ),	# radius
 	( 'Atm (km)'     , 9 ),    # atm thickness
	( 'R_s (km)'     , 10 ),    # safe radius (radius + max elevation or atm thickness)
	( 'mu (km^3/s^2)' , 11 ),	# gravitational parameter
	( 'mu_c (km^3/s^2)', 12 ),  # central body grav parameter
	( 'R_soi (km)'    ,  13 ), # SOI radius
	] )

   
body_names = OrderedDict( [
	( 0. , 'Moho'  ),
	( 1. , 'Eve'   ),
	( 2. , 'Kerbin'),
	( 3. , 'Duna'  ),      
	( 4. , 'Edna'  ),  
	( 5. , 'Dres'  ),      
	( 6. , 'Jool'  ),      
	( 7. , 'Lindor'),  
	( 8. , 'Eeloo' ),  
	( 9. , 'Hamek' ),  
	( 10., 'Nara'  ),  
	] )
	
# Defining NAMED TUPLES for accessing orbital elements by name
# >>> ['body'] + [ oe.split(' ')[0] for oe in col.keys()[1:] ]
# ['body', 'Epoch', 'a', 'e', 'i', 'Node', 'w', 'M', 'R', 'mu']
#


orbital_elements = namedtuple('orbital_elements', ['body'] + [ list(oe.split(' '))[0] for oe in list(col.keys())[1:] ] )
#orbital_elements = namedtuple('orbital_elements', ['body'] + [ oe.split(' ')[0] for oe in col.keys()[1:] ] )


body_tuple = OrderedDict( [
	( _body_name, orbital_elements( _body_name, *body[_brow][1:] ) )
	for (_body_name,_brow) in row.items()
	] )
	
#print(body_tuple['Moho'])

def _period( self ):
	"""
	Return's the body's orbital period, in days.
	"""
	return 2 * pi * sqrt( self.orbital_elements[0]**3 / self.mu_central_body )
	# http://en.wikipedia.org/wiki/Orbital_period#Calculation




# Instantiating PyKEP planet objects for each of the bodies (and making them indexable by body name)
# http://keptoolbox.sourceforge.net/documentation.html#PyKEP.planet

_planet = pk.planet.keplerian
### Problem is with the planet consutrctor or smth.

body_obj = OrderedDict( [
	( _b.body, _planet(
		# when
		# 	a PyKEP.epoch indicating the orbital elements epoch
		pk.epoch( _b.Epoch),#, pk.epoch.epoch_type.MJD ),
		
		# orbital_elements
		# 	a sequence of six containing a,e,i,W,w,M (SI units, i.e. meters and radiants)
		(	_b.a * 1000.,
			_b.e,
			radians( _b.i    ),
			radians( _b.Node ),
			radians( _b.w    ),
			_b.M,
			),
		
		# mu_central_body
		# 	gravity parameter of the central body (SI units, i.e. m^2/s^3)
		_b.mu_c * 1000.**3,
		# pk.MU_SUN == 1.32712428e+20 m^2/s^3
		
		# mu_self
		# 	gravity parameter of the planet (SI units, i.e. m^2/s^3)
		_b.mu * 1000.**3,	 # converting units: km^3/s^2 --> m^3/s^2
		
		# radius
		# 	body radius (SI units, i.e. meters)
		_b.R * 1000.,
		
		# safe_radius
		# 	body distance safe for a spacecraft fly-by
		_b.R_s * 1000.,
		# "The spacecraft cannot go below R_s at any time"
		
		# name
		# 	body name
		_b.body,
		) )
	for _b in body_tuple.values()
	] )


( Moho,Eve,Kerbin,Duna,Edna,Dres,Jool,Lindor,Eeloo,Hamek,Nara ) = body_obj.values()



# http://en.wikipedia.org/wiki/Orbital_elements
body_label = {
#	'Object UID'    : ,
	
	# a specified point in time
#	'Epoch (MJD)'   : ,
	
	# measure of the radius of an orbit taken from the points of that same orbit's two most distant points
	# http://en.wikipedia.org/wiki/Semimajor_axis
	'a (km)'        : r'$a$: semi major axis (km)',
	
	# The orbital eccentricity of an astronomical body is the amount by which its orbit deviates from a perfect circle, where 0 is perfectly circular, and 1.0 is a parabola, and no longer a closed orbit.
	# http://en.wikipedia.org/wiki/Orbital_eccentricity
	'e'             : r'$e$: eccentricity',
	
	# vertical tilt of the ellipse with respect to the reference plane, measured at the ascending node (where the orbit passes upward through the reference plane)
	# http://en.wikipedia.org/wiki/Inclination
	'i (deg)'       : r'$i$: inclination (deg)',
	
	# angle from a reference direction, called the origin of longitude, to the direction of the ascending node, measured in a reference plane
	# http://en.wikipedia.org/wiki/Longitude_of_the_ascending_node
	'Node (deg)'    : r'$\Omega$: longitude of the ascending node (deg)',
	
	# angle between the orbit's periapsis (the point of closest approach to the central point) and the orbit's ascending node (the point where the body crosses the plane of reference from South to North)
	# http://en.wikipedia.org/wiki/Argument_of_periapsis
	'w (deg)'       : r'$\omega$: argument of periapsis (deg)',
	
	# position of the orbiting body along the ellipse at a specific time (the "epoch")
	# http://en.wikipedia.org/wiki/Mean_anomaly
	'M (deg)'       : r'$M_o$: mean anomaly at epoch (deg)',
	
	# the body's radius	
	'R (km)'        : r'$R$: radius (km)',
	
	# product of the gravitational constant G and the mass M of the body
	# http://en.wikipedia.org/wiki/Standard_gravitational_parameter
	'mu (km^3/s^2)' : r'$\mu$: gravitational parameter (km^3/s^2)',
	}
	


#--------------------------------------# Problem constants

# 	Table 6: Other constants and conversions
#SECSINKDAY  = 12 * 60 * 60 
#DAYSINKYEAR = 365*12 * 60 * 60
