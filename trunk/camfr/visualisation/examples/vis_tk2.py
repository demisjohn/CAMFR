#!/usr/bin/env python

####################################################################
#
# Illustrates plotting from Tk.
#
####################################################################

from camfr import *
from camfr_tk import *
from Numeric import *

set_lambda(1)
set_N(50)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define slabs.

PML = -0.1

space = Slab(air(4 + 2*PML*1j))

# Loop over width.

v = [] # To keep track of the data points.

for W in arange(0.100, 0.200, 0.005):
    slab = Slab(air(2 + PML*1j - W/2.0) + GaAs(W)  \
              + air(2 + PML*1j - W/2.0))
    stack = Stack(space(0) + slab(0.5) + space(0))
    stack.calc()

    print W

    v.append((W,abs(stack.R12(1,1))))

    free_tmps()

plot_vector(v)    

