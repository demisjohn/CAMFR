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
set_precision(1000)
set_precision_enhancement(500)
set_dx_enhanced(.0001)

# Set geometry parameters

GaAs = Material(3.4)
wg_m = Material(1.8)
air  = Material(1.0)
  
a = .600     # period
r = .150/2.0 # rod radius

set_left_wall(slab_H_wall)

PML = 0

cl = 0 # air cladding

periods = 3  # periods above outer waveguide
sections = 1 # intermediate 90 deg sections

# Define slabs.

no_rods = Slab(air(a-r+(sections+1+periods)*a+cl-PML*1j))

# Central waveguide.
 
cen = Slab(  air(a-r)                                               \
           + (sections+1+periods)*(GaAs(2*r) + air(a-2*r))          \
           + air(cl-PML*1j) )

# Vertical section.

ver = Slab(  air(a-r + (sections+1)*a)                              \
           + periods*(GaAs(2*r) + air(a-2*r) )                      \
           + air(cl-PML*1j) )

# Outer arms.
 
arm = Slab(  GaAs(r) + air(a-2*r)                                   \
           + sections*(GaAs(2*r) + air(a-2*r))                      \
	   + air(a)                                                 \
           + periods*(GaAs(2*r) + air(a-2*r))                       \
           + air(cl-PML*1j) )

# Find lowest order waveguide mode.

wg = BlochStack(cen(2*r) + no_rods(a-2*r))
wg.calc()

print wg

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

for x in arange(0.000, no_rods.width().real - cl - a, a/20.):
    for z in arange(0.000, splitter.length().real, a/20.):
        print >> outfile, abs(splitter.field(Coord(x, 0, z)).E2()),
    print >> outfile

outfile.close()    

