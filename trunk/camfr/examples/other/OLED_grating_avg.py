#!/usr/bin/env python

######################################################################
#
# Extraction efficiency in an organic LED with gratings
# averaged over dipole position.
#
######################################################################

from GARCLED import *

def calc(L, h, D, orders, wavelength=0.565):

  # Set parameters.

  set_period(L, L)
  
  set_lambda(wavelength)

  set_fourier_orders(orders, orders)

  print
  print "wavelength, orders L h:", wavelength, orders, L, h

  # Create materials.
  
  Al    = Material(1.031-6.861j)
  Alq3  = Material(1.655)
  NPD   = Material(1.807)
  ITO   = Material(1.806-0.012j)
  SiO2  = Material(1.48)
  SiN   = Material(1.95)  
  glass = Material(1.528)
  air   = Material(1)

  eta_sub, eta_out = [], []

  # Define layer structure.

  for pos in ['1', '2', '3']:

    print
    print "*** Position", pos, "***"
    print
  
    top = Uniform(Alq3,  0.050) + \
          Uniform(Al,    0.150) + \
          Uniform(air,   0.000)
  
    bot = Uniform(Alq3,  0.010) + \
          Uniform(NPD,   0.045) + \
          Uniform(ITO,   0.050) + \
          SquareGrating(D, SiO2, SiN, pos, h) + \
          Uniform(glass, 0.000)
  
    sub = Uniform(glass, 0.000) + \
          Uniform(air,   0.000)

    ref = Uniform(3) # Don't choose too high: otherwise takes too long
  
    cav = RCLED(top, bot, sub, ref)

    # Calculate average over source orientation for a given position.
     
    res = cav.calc(sources = [vertical, horizontal_x, horizontal_y],
                   steps=100, symmetric=True)

    eta_sub.append(res.eta_sub)
    eta_out.append(res.eta_out)
      
    sys.stdout.flush()
  
    free_tmps()

  # Print average over position and source.

  print
  print "*** Averaged over dipole position and source orientation ***"
  print
  print "Averaged extraction efficiency to substrate            :", \
                                                             average(eta_sub) 
  print "Averaged extraction efficiency to bottom outside world :", \
                                                             average(eta_out)

# Loop over period and thickness of the grating.

for L in [.300, .500, .700]:
  for h in [.100, .200, .400]:
    calc(L, h, L/2, orders=3)
  
