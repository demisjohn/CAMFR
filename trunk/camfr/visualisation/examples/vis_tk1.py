#!/usr/bin/env python

####################################################################
#
# Simple visualisation example.
#
####################################################################

from camfr import *

set_lambda(1)
set_N(50)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define waveguide.

set_lower_PML(-0.1)
set_upper_PML(-0.1)

slab = Slab(air(2) + GaAs(0.5) + air(2))

slab.calc()

# Visualise.

r_x = arange(1.0, 3.5, 0.01)

print "Plotting mode distribution (close window to continue)"
plot_neff(slab)

print "Plotting E field (close window to continue)"
plot_field(slab.mode(0), lambda f : f.E2().real, r_x)

print "Plotting H field (close window to continue)"
plot_field(slab.mode(0), lambda f : f.H2().real, r_x)

print "Plotting refractive index profile (close window to continue)"
plot_n(slab, r_x)
