#!/usr/bin/env python

#############################################################################
#
# Simple waveguide plotting example
#
#############################################################################

from camfr import *

set_lambda(1)

set_N(20)

set_lower_PML(-0.5)
set_upper_PML(-0.7)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define waveguide.

slab = Slab(air(2) + GaAs(0.5) + air(2))

slab.calc()

# Plot out some waveguide characteristics. 

plot(slab)

