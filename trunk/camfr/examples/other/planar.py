#!/usr/bin/env python

####################################################################
#
# planar 1D DBR
#
####################################################################

from camfr import *
from Numeric import *

# Define materials.

set_lambda(1.)
GaAs = Material(3.53)
AlAs = Material(2.95)

d_GaAs = .25*GaAs.n().real
d_AlAs = .25*AlAs.n().real

# TE incidence.

set_polarisation(TE)

GaAs_TE = Planar(GaAs)
AlAs_TE = Planar(AlAs)

s_TE = Stack(GaAs_TE(0) + 10*(GaAs_TE(d_GaAs) + AlAs_TE(d_AlAs)) + GaAs_TE(0))
GaAs_TE.set_theta(0)

for l in arange(0.8, 1.2, 0.05):
  set_lambda(l)
  s_TE.calc()
  print l, abs(s_TE.R12(0,0))**2
