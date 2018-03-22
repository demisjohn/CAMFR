#!/usr/bin/env python

####################################################################
#
# Parabolic refractive index profile.
#
####################################################################

from camfr import *

set_lambda(1)
set_N(20)
set_polarisation(TE)

w = 5.0  # width of waveguide

set_upper_PML(-0.1)
set_lower_PML(-0.1)

# Define parabolic refractive index profile.

def index(x):
    max_n = 3.5  # max refractive index
    a = 0.1      # slope
    n = max_n - a * pow(w / 2.0 - x, 2)
    if n < 1:
	return 1
    else:
	return n

# Construct a staircase approximation.

expr = Expression()
materials = [] 
steps = 10
for i in range(steps):
    x = i * w / steps
    m = Material(index(x + 0.5 * w / steps))
    materials.append(m)
    d = w / steps
    expr.add(m(d))

slab = Slab(expr)

# Compare continuous and staircase profile.

outfile = file("tutorial7.out",'w')

steps2 = 100
for i in range(steps2):
    x = i * w / steps2
    #print >> outfile, x, index(x), slab.n(Coord(x, 0, 0)).real
        outfile.writelines(  [str(x), str(index(x)), str(slab.n(Coord(x, 0, 0)).real) ]  )


outfile.close()
