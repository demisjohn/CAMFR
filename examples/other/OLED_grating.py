#!/usr/bin/env python

####################################################################
#
# Extraction efficiency in an organic LED with gratings
#
####################################################################

from GARCLED import *

# Set parameters.
  
set_lambda(0.565)

orders = 2
set_fourier_orders(orders,orders)

L = 0.400
set_period(L, L)

# Create materials.
  
Al    = Material(1.031-6.861j)
Alq3  = Material(1.655)
NPD   = Material(1.807)
SiO2  = Material(1.48)
SiN   = Material(1.95) 
ITO   = Material(1.806-0.012j)
glass = Material(1.528)
air   = Material(1)

# Define layer structure.
 
top = Uniform(Alq3,  0.050) + \
      Uniform(Al,    0.150) + \
      Uniform(air,   0.000)
  
bot = Uniform(Alq3,  0.010) + \
      Uniform(NPD,   0.045) + \
      Uniform(ITO,   0.050) + \
      SquareGrating(L/4, SiO2, SiN, "1", 0.100) + \
      Uniform(glass, 0.000)
  
sub = Uniform(glass, 0.000) + \
      Uniform(air,   0.000)

ref = Uniform(3) # Don't choose too high: otherwise takes too long

cav = RCLED(top, bot, sub, ref)

# Calculate.
  
cav.calc(sources=[vertical, horizontal_x], weights=[1, 2],
         steps=30, symmetric=True)

cav.radiation_profile(horizontal_x, steps=30)
