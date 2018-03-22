#!/usr/bin/env python

####################################################################
#
# Another stack example
#
####################################################################

from camfr import *

# Set constants.

set_lambda(1.5)
set_N(20)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define waveguide sections.

set_lower_PML(-0.1)
set_upper_PML(-0.1)

normal = Slab(air(2.0) + GaAs(0.5) + air(2.0))
thick  = Slab(air(1.9) + GaAs(0.7) + air(1.9))

# Calculate reflection of the fundamental mode for different 
# lengths of the central thick section.

outfile = file("tutorial3.out",'w')

for L in arange(0.000, 0.500, 0.010):    
    stack = Stack(normal(0) + thick(L) + normal(0))
    stack.calc()
    #print >> outfile, L, abs(stack.R12(0,0))
    outfile.writelines(  [ str( L ), str( abs(stack.R12(0,0)) ) ]  )


outfile.close()
