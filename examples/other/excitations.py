#! /usr/bin/env python

###################################################################
#
# Illustrates different excitations for cartesian stacks.
#
###################################################################

from camfr import *

set_N(50)
set_lambda(1.55)

GaAs = Material(3.5)
air  = Material(1.0)

slab = Slab(air(2-.1j)+ GaAs(1)+ air(2-.1j))
s = Stack(slab(1))

eps = 1e-3 # Precision for calculating overlap integrals.

# General excitation using a Python function.

A     = 1.0
x0    = slab.width().real/2.
sigma = 0.5

def f(x):
  return A*exp(-0.5*((x-x0)/sigma)**2)

s.set_inc_field_function(f, eps)

# Faster variant using a built-in Gaussian excitation.

s.set_inc_field_gaussian(A, sigma, x0, eps)

# Plane wave: f = A + slope*x

slope = 0.0
s.set_inc_field_plane_wave(A, slope, eps)
