#!/usr/bin/env python

# Author: Mihai Ibanescu <michel@mit.edu>
#
# This example file simulates transmission between two identical omniguides
# that have an air gap separating them.
#
# Transmission from the TE01 mode of the first waveguide into the TE01 mode
# of the second waveguide is calculated for various values of zsep, the
# length of the air gap.
# 
# A field animation is also included.

from camfr import *

set_lambda(1/0.35)
set_N(20)

set_circ_order(0)
set_polarisation(TE)
set_circ_PML(-0.2)

air  = Material(1.0)
mat1 = Material(3.0) 
mat2 = Material(1.5)

Rcore = 4.0 # core radius
Nbi = 2     # number of bilayers
a = 1.0     # bilayer thickness
d1 = 0.30   # mat1 thickness
d2 = a - d1 # mat2 thickness
p = 2.0;    # air padding

# Define omniguide and gap.

og = Circ(air(Rcore) + Nbi*(mat1(d1)+mat2(d2)) + mat1(d1) + air(p))
space = Circ(air(Rcore+Nbi*a+d1) + air(p))

og.calc()
te01 = 0
print "Omniguide modes:"
for i in arange(0, N(), 1):
    print og.mode(i).n_eff()
    # Find TE01 mode as having n_eff closest to 1.
    if abs(1 - og.mode(i).n_eff().real) < abs(1 - og.mode(te01).n_eff().real):
        te01 = i

# Excite TE01 mode only.

print " zsep   Transmission"
for zsep in arange(0,4.01,0.2):
    s = Stack(og(0.0) + space(zsep) + og(0.0))
    s.calc()
    print zsep, abs(s.T12(te01,te01))**2

inc = zeros(N())
inc[te01] = 1
zsep = 2.0
s = Stack(og(0.0) + space(zsep) + og(0.0))
s.set_inc_field(inc)
s.calc()

vx = arange(0,Rcore+Nbi*a+d1+p,0.15)
vz = arange(-3,8,0.15)

animate_field(s, lambda f : f.E2(), vx, vz)
