#! /usr/bin/env python

###################################################################
#
# Calculates modification of spontaneous emission of a dipole
# between two metal plates.
#
###################################################################

from camfr import *

set_N(60)
set_lambda(1)
set_circ_order(1)
set_circ_field_type(cos_type)

# Define waveguide and wall.

set_circ_PML(-0.5)

air_m = Material(1.0)
air = Circ(air_m(10))
air.calc()

wall = E_Wall(air)

# Calculate change in spontaneous emission rate.

for d in arange(0.01, 3.0, 0.05):     

    # Define cavities.

    half      = Stack(air(d/2.) + wall)
    half_open = Stack(air(d/2.))

    source_pos  = Coord(0,0,0)
    orientation = Coord(1,0,0)
    
    cav = Cavity(half, half)
    cav.set_source(source_pos, orientation)

    cav_open = Cavity(half_open, half_open)
    cav_open.set_source(source_pos, orientation)

    # Analytic formula for spontaneous emission rate.

    x = floor(d+0.5)
    exact = 3.*x/4./d + pow(x/d,3)/4. - x/16./d/d/d

    # Numerical formula as ratio of total emitted powers.

    numeric =   half.     field(Coord(0,0,0)).E1().real   \
	      / half_open.field(Coord(0,0,0)).E1().real 

    print d, exact, numeric

