#! /usr/bin/env python

##############################################################################
#
# Simple BlochStack plotting example.
#
##############################################################################

from camfr_work import *

set_lambda(1.5)
set_N(60)

# Set geometry parameters.

GaAs = Material(3.4)
air  = Material(1.0)
  
a = .600     # period
r = .150/2.0 # rod radius
periods = 5  # periods in lateral direction

set_lower_wall(slab_H_wall)

# Define slabs.

no_rods = Slab(air(a-r+periods*a))
rods    = Slab(air(a-r) + periods*(GaAs(2*r) + air(a-2*r)) )

# Plot BlochStack.

wg = BlochStack(rods(2*r) + no_rods(a-2*r) + rods(2*r) + no_rods(a-2*r))
wg.calc()

plot(wg)

