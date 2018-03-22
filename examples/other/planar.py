#!/usr/bin/env python

####################################################################
#
# planar 1D DBR
#
####################################################################

from camfr import *

# Define structure.

set_lambda(1)

GaAs_m = Material(3.5)
AlAs_m = Material(2.9)

GaAs = Planar(GaAs_m)
AlAs = Planar(AlAs_m)

d_GaAs = get_lambda()/4./GaAs_m.n().real
d_AlAs = get_lambda()/4./AlAs_m.n().real

s = Stack(GaAs(0) + 10*(GaAs(d_GaAs) + AlAs(d_AlAs)) + GaAs(0))

# Loop over incidence angles.

for theta in arange(0, 90, 0.5):

  GaAs.set_theta(theta * pi / 180.)
  print(theta)
  
  set_polarisation(TE)
  s.calc()
  print(abs(s.R12(0,0))**2 )

  set_polarisation(TM)
  s.calc()
  print(abs(s.R12(0,0))**2)
  
