#! /usr/bin/env python

###################################################################
#
# A semi-infinite photonic crystal.
#
###################################################################

from camfr import *
from Numeric import *

set_lambda(1.5)
set_N(50)

# Set geometry parameters

GaAs = Material(3.4)
air  = Material(1.0)
  
a = .600     # period
r = .150/2.0 # rod radius

set_lower_wall(slab_H_wall)

cl = 0       # air cladding
periods = 4  # lateral periods

# Define slabs.

inc_wg = Slab(GaAs(1.5*r) + air(a-2.5*r+periods*a+cl))

no_rods = Slab(air(a-r+periods*a+cl))
 
cen = Slab(  air(a-r)                                               \
           + periods*(GaAs(2*r) + air(a-2*r))                       \
           + air(cl) )

# Calculate semi-infinite stack.

s_inf = InfStack(cen(2*r) + no_rods(a-2*r))
s = Stack(inc_wg(a) + s_inf)

s.calc()

print abs(s.R12(0,0))**2
