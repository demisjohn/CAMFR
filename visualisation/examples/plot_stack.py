#! /usr/bin/env python

#############################################################################
#
# Simple stack plotting example.
#
#############################################################################

from camfr import *

set_polarisation(TM)
set_N(100)
set_lambda(1.55)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define geometry.

wg  = Slab(air(0.5) + GaAs(0.2) + air(1.0))
gap = Slab(air(1.7))

s = Stack(wg(1) + gap(1) + wg(2))

# Set incident field.

inc = zeros(N())
inc[0] = 1
s.set_inc_field(inc)

# Do some plotting.

plot(s)

