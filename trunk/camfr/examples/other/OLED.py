#!/usr/bin/env python

####################################################################
#
# Extraction efficiency in an organic LED
#
####################################################################

from RCLED import *

# Set parameters.

set_lambda(0.565)

# Create materials.
  
Al    = Material(1.031-6.861j)
Alq3  = Material(1.655)
NPD   = Material(1.807)
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
      Uniform(glass, 0.000)

sub = Uniform(glass, 0.000) + \
      Uniform(air,   0.000)

ref = Uniform(10)

cav = RCLED(top, bot, sub, ref)

# Calculate.

cav.calc()

cav.radiation_profile(horizontal)

