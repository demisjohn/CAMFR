#!/usr/bin/env python

####################################################################
#
# Simple waveguide example
#
####################################################################

from camfr import *

set_lambda(1)
set_N(20)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define waveguide.

slab = Slab(air(2) + GaAs(0.5) + air(2))

slab.calc()

# Print out some waveguide characteristics. 

print( slab.mode(0).kz() )
print( slab.mode(1).n_eff() )
print( slab.mode(2).field(Coord(2.25, 0, 0)) )
print( slab.mode(3).field(Coord(2.25, 0, 0)).E2() )
print( slab.mode(4).field(Coord(2.25, 0, 0)).E2().real )

# Do some interactive plotting.

slab.plot()
