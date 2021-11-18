# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:56:18 2019

@author: v.fave
"""

"""
Defines Kerbol system in terms of pykep planets.

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


# Values from KSP wiki
# 	Table 4: Keplerian orbit elements of the kerbol planets at t=0
# 	Table 5: Satellite physical constants
# 	Table 6: Other constants and conversions
G = 6.674e-11 #* 1.000027

body = np.array( [
	# UID   Epoch    a (km)             e           i (deg)       Node (deg)    w (deg)     M (rad)  R (km)      Atm_h     Rs (km)    mu (km^3/s^2)     mu_c              R SOI (km)   
	[  0.,  0.,   5263138.304 ,        0.2,          7.0,           70.0,        15.0,      3.14,     250.0,      0.0,      257.0,     1.6860938e2,  	1.756567e+19*G,    9646.663     ],
	[  1.,  0.,   9832684.544 ,        0.01,         2.1,           15.0,         0.0,      3.14,     700.0,     90.0,      791.0,     8.1717302e3,  	1.756567e+19*G,    85109.365    ],
	[  2.,  0.,       31500.0 ,        0.55,        12.0,           80.0,        10.0,      0.90,      13.0,      0.0,       20.0,     8.2894498e-3,   8.1717302e3,        126.12327    ],
	[  3.,  0.,  13599840.256 ,        0.0,          0.0,            0.0,         0.0,      3.14,     600.0,     70.0,      671.0,     5.2915793e13*G, 1.756567e+19*G,     84159.286    ],    
	[  4.,  0.,     12000.0   ,        0.0,          0.0,            0.0,         0.0,      1.40,     200.0,      0.0,      208.0,     6.5138398e1,    3.5316000e3,        2429.5591    ], 
	[  5.,  0.,     47000.0   ,        0.0,          6.0,           38.0,        78.0,      0.90,      60.0,      0.0,       66.0,       1.7658000,    3.5316000e3,        2247.4284    ],  
	[  6.,  0.,  20726155.264 ,        0.051,        0.06,         135.5,         0.0,      3.14,     320.0,     50.0,      371.0,     3.0136321e2,  	1.756567e+19*G,     47921.949   ],       
	[  7.,  0.,      3200.000 ,        0.03,         0.2,            0.0,         0.0,       1.7,     130.0,      0.0,      143.0,     1.8568369e1,    3.0136321e2,        1049.5989    ], 
	[  8.,  0.,  40839348.203 ,        0.145,        5.0,          280.0,        90.0,      3.14,     138.0,      0.0,      143.0,     2.1484489e1,  	1.756567e+19*G,     32832.840   ],     
	[  9.,  0.,  68773560.320 ,        0.05,         1.304,          52.0,        0.0,      0.10,    6000.0,    200.0,     6200.0,     2.8252800e5,  	1.756567e+19*G,     2.4559852e9 ],     
	[  10., 0.,  90118820.000 ,        0.26,         6.15,          50.0,       260.0,      3.14,     210.0,      0.0,      215.0,     1.1149358e2,    1.756567e+19*G,     1.1908294e5  ],      
	], dtype=np.float64 )
	 

#Kerbin = pk.planet.keplerian(pk.epoch(0), (13599840256 ,0,0,0,0,3.14), 1.1723328e18, 3.5316000e12, 600000, 670000 , 'Kerbin')
#Duna   = pk.planet.keplerian(pk.epoch(0), (20726155264 ,0.051,deg2rad(0.06) ,0,deg2rad(135.5),3.14), 1.1723328e18, 3.0136321e11,320000, 370000 , 'Duna')
# Mapping from name of a body to number of row in which it is represented in the body matrix
row = OrderedDict( [
	( 'Moho'   , 0 ),
	( 'Eve'    , 1 ),
	( 'Gilly'  , 2 ),
	( 'Kerbin' , 3 ),
	( 'Mun'    , 4 ),
	( 'Minmus' , 5 ),
	( 'Duna'   , 6 ),
	( 'Ike'    , 7 ),
	( 'Dres'   , 8 ),    
	( 'Jool'   , 9 ),    
	( 'Eeloo'  , 10 ),    
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
	( 1. , 'Eve'  ),
	( 2. , 'Gilly'  ),
	( 3. , 'Kerbin'  ),
	( 4. , 'Mun'  ),
	( 5. , 'Minmus'  ),
	( 6. , 'Duna'  ),    
	( 7. , 'Ike'  ),      
	( 8. , 'Dres'  ),      
	( 9. , 'Jool'  ),        
   ( 10., 'Eeloo'  ),         
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


( Moho,Eve,Gilly,Kerbin,Mun,Minmus,Duna,Ike,Dres, Jool, Eeloo ) = body_obj.values()



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
SECSINKDAY  = 6 * 60 * 60 
DAYSINKYEAR = 2556.5 * 60 * 60
