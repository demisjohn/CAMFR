#!/usr/bin/env python

####################################################################
#
# Simple example of the visualisation using Matlab.
#
####################################################################

from camfr import *
from camfr_matlab import *

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

print "Plotting mode distribution"
plot_neff(slab)
raw_input("Press <enter> to continue")

print "Plotting E field"
plot_field(slab, modenumber, lambda f : f.E2().real)
raw_input("Press <enter> to continue")

print "Plotting H field"
plot_field(slab, modenumber, lambda f : f.H2().real)
raw_input("Press <enter> to continue")

print "Plotting refractive index profile"
plot_n(slab)
raw_input("Press <enter> to continue")
