#!/usr/bin/env python

####################################################################
#
# Calculate field profiles in a stack, but exploit symmetry
#
####################################################################

from camfr import *
from Numeric import *

set_lambda(1)
set_N(10)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define stack.

PML = 0.1

set_left_wall(slab_H_wall)

slab  = Slab(GaAs(.25) + air(2 - PML*1j))
space = Slab(air(2.25 - PML*1j))

stack = Stack(space(0) + slab(0.5) + space(0))

# Set incident field and calculate stack.

inc = zeros(N())
inc[0] = 1
stack.set_inc_field(inc)

# Save the field to a file.

outfile = open("tutorial6.out",'w')

for x in arange(0.000, 2.250, 0.100):
    for z in arange(0.000, 0.500, 0.010):
        print >> outfile, abs(stack.field(Coord(x, 0, z)).E2()),
    print >> outfile

outfile.close()

