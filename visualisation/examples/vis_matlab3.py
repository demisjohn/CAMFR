#!/usr/bin/env python

####################################################################
#
# Illustrates passing data to and from Matlab.
#
####################################################################

from camfr import *
from Numeric import *
from camfr_matlab import *

set_lambda(1)
set_N(50)
set_polarisation(TE)

# Define materials.

GaAs = Material(3.5)
air  = Material(1.0)

# Define stack.

PML = 0.1

space = Slab(air(4 - 2*PML*1j))
slab = Slab(air(1.5 - PML*1j) + GaAs(1) + air(1.5 - PML*1j))

stack = Stack(space(0) + slab(0.5) + space(0))

stack.calc()

# Use Matlab to do some manipulations on a matrix (Note: this can 
# also be done directly using Numerical Python).

matlab.put("R12", stack.R12())
matlab("A = R12 * R12")
R12_R12 = matlab.get("A")

# Illustrates multiline commands using triple quotes

matlab("""
surf(abs(R12))
figure
pcolor(abs(R12*R12))
""")

raw_input("Press <enter> to continue")
