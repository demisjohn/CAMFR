#!/usr/bin/env python

####################################################################
#
# Simple stack example
#
####################################################################

from camfr import *
from Numeric import *

set_lambda(1)
set_N(20)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define slabs.

slab  = Slab(air(2) + GaAs(0.5) + air(2))
space = Slab(air(4.5))

# Print the reflectivity for different lengths.

for L in arange(0.005, 0.100, 0.005):
    stack = Stack(space(0) + slab(L) + space(0))
    stack.calc()
    print L, abs(stack.R12(0,0))
