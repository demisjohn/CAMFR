#!/usr/bin/env python

####################################################################
#
# Calculate field profiles in a stack.
#
####################################################################

from camfr import *
from Numeric import *

set_lambda(1)
set_N(20)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define stack.

PML = 0.1

slab  = Slab(air(2 - PML*1j) + GaAs(.5) + air(2 - PML*1j))
space = Slab(air(4.5 - 2*PML*1j))

stack = Stack(space(0) + slab(0.5) + space(0))

# Set incident field and calculate stack.

inc = zeros(N())
inc[0] = 1
stack.set_inc_field(inc)

stack.calc()

# Save the field to a file.

outfile = open("tutorial5.out",'w')

for x in arange(0.000, 4.500, 0.100):
    for z in arange(0.000, 0.500, 0.010):
	print >> outfile, x, z, abs(stack.field(Coord(x - PML*1j, 0, z)).E2())
    print >> outfile

outfile.close()
