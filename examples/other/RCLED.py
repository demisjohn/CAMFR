#!/usr/bin/env python

####################################################################
#
# Calculates spontaneous emission rate and extraction efficiency in
# planar RCLEDs.
#
####################################################################

from camfr_work import *
from cmath import *
from Numeric import *

set_lambda(0.980)
set_circ_order(1)

# Define materials.

metal_m = Material(0.2-6.5j)
air_m   = Material(1.0)
GaAs_m  = Material(3.5)
QW_m    = Material(3.5-0.1j)
AlOx_m  = Material(1.55)

d_AlOx = get_lambda()/4.0/AlOx_m.n().real
d_GaAs = get_lambda()/4.0/GaAs_m.n().real

# Define the calculate function

def calc(r, PML, M, max_rho, spacer, resolution):
  
  # 1-D structures

  set_N(M)

  set_circ_PML(PML)
  
  metal = Circ(metal_m(r))
  air   = Circ(  air_m(r)) 
  GaAs  = Circ( GaAs_m(r))
  QW    = Circ(   QW_m(r))
  AlOx  = Circ( AlOx_m(r))

  # Reference bulk emitter

  half_open = Stack(GaAs(1))

  r0 = Coord(0,0,0)
  orientation = Coord(1,0,0)

  open = Cavity(half_open,half_open)
  open.set_source(r0, orientation)

  ref_rate = (half_open.field(r0).E1()).real
  
  # 2-D structure

  top = Stack(GaAs(0.0)+QW(.0025)+GaAs(.048-.0025)+metal(.120)+air(.002))
  bot = Stack(GaAs(0.0)+QW(.0025)+GaAs(spacer)+AlOx(d_AlOx)+GaAs(0.0))
  
  # Cavity

  rcled = Cavity(bot,top)
  rcled.set_source(r0, orientation)

  # Calculate spontaneous emission rate.

  rate = rcled.field(r0).E1().real / ref_rate

  # Calculate QW emission.

  S_top_QW = top.inc_S_flux(0, max_rho, resolution)
  S_bot_QW = bot.inc_S_flux(0, max_rho, resolution)

  # Calculate bottom exit emission.
  # Assume incoherent reflection from bottom substrate-air interface.
  
  outcoupling = Stack(GaAs(0.0) + air(0.0))
  outcoupling.set_inc_field(bot.trans_field())

  S_top_out =         top.ext_S_flux(0, max_rho, resolution)
  S_bot_out = outcoupling.ext_S_flux(0, max_rho, resolution)

  # Calculate extraction efficiency.

  eta = S_bot_out / (S_bot_QW + S_top_QW)

  # Crude approximate of coupling efficiency to NA.

  NA = 0.5
  out = outcoupling.trans_field()
  k = 2.0 * pi / get_lambda() # n = 1 in air

  P = P_NA = 0
  for i in range(N()):
    kz = air.mode(i).kz()
    theta = acos(kz.real / k).real

    if abs(kz.real) > abs(kz.imag):
      P_theta = pow(abs(out[i]), 2) * sin(theta)
      P += P_theta

      if (sin(theta) <= NA):
        P_NA += P_theta
    
  # Report results

  print "@ ", spacer, r, PML, M, max_rho, resolution, ":", \
        rate, eta, P_NA/P, eta*P_NA/P, \
        S_top_QW, S_bot_QW, S_top_out, S_bot_out

  free_tmps()



# Main loop.

eps = 1e-10

for spacer in arange(0.130, 0.160, 0.002):
  for r in arange(10.0, 10.0+eps, 2.0):
    for PML in arange(-0.05, -0.05+eps, 0.01):
      for M in arange(300, 301, 20):
        max_rho = r
        resolution = 1e-10
        calc(r, PML, M, max_rho, spacer, resolution)
        
    
