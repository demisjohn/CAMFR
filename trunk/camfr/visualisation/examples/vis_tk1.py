#!/usr/bin/env python

####################################################################
#
# Simple example of the visualisation based on the Tk toolkit.
#
####################################################################

from camfr import *
from camfr_tk import *

set_lambda(1)
set_N(50)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define waveguide.

PML = -0.1
slab = Slab(air(2+PML*1j) + GaAs(0.5) + air(2+PML*1j))

slab.calc()

# Visualise.

modenumber = 0

print "Plotting mode distribution (close window to continue)"
plot_neff(slab)

print "Plotting E field (close window to continue)"
plot_field(slab, modenumber, lambda f : f.E2().real)

print "Plotting H field (close window to continue)"
plot_field(slab, modenumber, lambda f : f.H2().real)

print "Plotting refractive index profile (close window to continue)"
plot_n(slab)
