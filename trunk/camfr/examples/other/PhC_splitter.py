#! /usr/bin/env python

###################################################################
#
# 90 deg 3dB splitter in rectangular lattice of rectangular 
# GaAs rods in air
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

cl = 0 # air cladding

periods = 3  # periods above outer waveguide
sections = 1 # intermediate 90 deg sections

# Define slabs.

no_rods = Slab(air(a-r+(sections+1+periods)*a+cl))

# Central waveguide.
 
cen = Slab(  air(a-r)                                               \
           + (sections+1+periods)*(GaAs(2*r) + air(a-2*r))          \
           + air(cl) )

# Vertical section.

ver = Slab(  air(a-r + (sections+1)*a)                              \
           + periods*(GaAs(2*r) + air(a-2*r) )                      \
           + air(cl) )

# Outer arms.
 
arm = Slab(  GaAs(r) + air(a-2*r)                                   \
           + sections*(GaAs(2*r) + air(a-2*r))                      \
	   + air(a)                                                 \
           + periods*(GaAs(2*r) + air(a-2*r))                       \
           + air(cl) )

# Find lowest order waveguide mode.

wg = BlochStack(cen(2*r) + no_rods(a-2*r))
wg.calc()

guided = 0
for i in range(2*N()):
    if (abs(wg.mode(i).kz().imag) < abs(wg.mode(guided).kz().imag)):
        if wg.mode(i).kz().real > 0:
            guided = i

# Calculate splitter.

splitter = Stack(  5*(cen(2*r) + no_rods(a-2*r))                   \
                 +    ver(2*r) + no_rods(a-2*r)                    \
                 + 5*(arm(2*r) + no_rods(a-2*r)) )

splitter.set_inc_field(wg.mode(guided).fw_field())
splitter.calc()

print "R", splitter.R12(0,0)

# Calculate field.

outfile = open("splitter.out", 'w')

for x in arange(0.000, no_rods.width() - cl - a, a/20.):
    for z in arange(0.000, splitter.length(), a/20.):
        print >> outfile, abs(splitter.field(Coord(x, 0, z)).E2()),
    print >> outfile

outfile.close()    

