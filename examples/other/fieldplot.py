#! /usr/bin/env python

####################################################################
#
# Illustrates 2D.
#
####################################################################

from camfr import *

set_N(40)
set_lambda(1.55)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define geometry.

set_lower_PML(-0.1)
set_upper_PML(-0.1)

wg  = Slab(air(2.0) + GaAs(0.2) + air(2.0))
gap = Slab(air(4.2))

s = Stack(wg(0) + gap(1) + wg(0))

inc = zeros(N())
inc[0] = 1
s.set_inc_field(inc)

s.calc()

# Do some plotting.

r_x = arange( 1.5, 2.8, 0.05)
r_z = arange(-1.0, 2.0, 0.05)

print "Doing some plots. Close window to continue to the next one."

wg.plot_n(r_x)

wg.mode(0).plot_field(lambda f : f.E2().real, r_x)

s.plot_field(lambda f : f.E2().real, r_x, r_z)

s.animate_field(lambda f : f.E2(), r_x, r_z)

