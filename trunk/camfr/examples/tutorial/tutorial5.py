#!/usr/bin/env python

####################################################################
#
# Calculate field profiles in a stack.
#
####################################################################

from camfr import *

set_lambda(1)
set_N(20)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define stack.

set_lower_PML(-0.1)
set_upper_PML(-0.1)

slab  = Slab(air(2) + GaAs(.5) + air(2))
space = Slab(air(4.5))

stack = Stack(space(0) + slab(0.5) + space(0))

# Set incident field and calculate stack.

inc = zeros(N())
inc[0] = 1
stack.set_inc_field(inc)

stack.calc()

# Do some interactive plotting.

stack.plot()
