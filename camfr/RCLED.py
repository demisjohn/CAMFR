#! /usr/bin/env python

from camfr import *

#############################################################################
#
# Simple library to simulate planar RCLEDs.
#
# Version 3.0: 20070705
#
#############################################################################

#############################################################################
#
# Pick correct sign for square root.
#
#############################################################################

def safe_sqrt(k2):

  k = sqrt(k2)
  
  if k.imag > 0: 
    k = -k

  if abs(k.imag) < 1e-9:
    if (k.real < 0):
      k = -k

  return k



#############################################################################
#
# RCLED class
#
# top and bot: Expressions describing the cavity, from the source outwards.
#
# sub: Expression describing the substrate placed under the bottom part,
# through which the power is extracted incoherently and in a single pass.
#
# ref: Planar describing an infinitely thin reference layer
# where the emission occurs. If the index of this layer is high enough, then
# all the plane waves emitted by the source will be propgating in the
# reference layer, which makes the numerics easier. There are however
# normalisation factors (powers of self.corr) to be taken into account
# for the total power then. See calc() function.
#
#############################################################################

class RCLED:

  def __init__(self, top, bot, sub, ref = None):

    if ref:
      self.cor = ref.core().n() / top.inc().core().n()
      self.top = Stack(ref(0) + top)
      self.bot = Stack(ref(0) + bot)
    else:
      self.cor = 1.0
      self.top = Stack(top)
      self.bot = Stack(bot)

    self.sub = Stack(sub)

    # Consistency check.

    if self.top.inc().core().n() != self.bot.inc().core().n():
      raise Exception("Incidence media of top and bottom stack don't match.")

    if self.bot.ext().core().n() != self.sub.inc().core().n():
      raise Exception("Exit medium of bottom doesn't match inc medium of sub.")

  # Todo: move the function code itself here to get better design.
  
  def calc(self, sources=None, weights=None):
    return calc(self, sources, weights)

  def radiation_profile(self, source, steps=250, location='out',
                        density='per_solid_angle', show=True, filename=None):
    return radiation_profile(self,source,steps,location,density,show,filename)



#############################################################################
#
# Source terms for vertical dipole, i.e. perpendicular to layer interfaces.
# These are given per unit solid angle.
#
#############################################################################

def vertical(kt, cav):

  if get_polarisation() == TM:

    n = cav.top.inc().core().n()
    k0 = 2.*pi/get_lambda()
    A  = sqrt(3./8./pi) * kt/k0/n
    return A, -A

  else:

    return 0.0, 0.0



#############################################################################
#
# Source terms for horizontal dipole with averaged azimuth.
# These are given per unit solid angle.
#
# Note: TE carries 75% of the total power, TM 25%.
#
#############################################################################

def horizontal(kt, cav):
  
  if get_polarisation() == TE:

    A = sqrt(3./16./pi)
    return A, A

  else:

    n = cav.top.inc().core().n()
    k0 = 2.*pi/get_lambda()
    kz = safe_sqrt(k0*k0*n*n - kt*kt)
    A  = sqrt(3./16./pi) * kz/k0/n
    return A, A



#############################################################################
#
# Calculate power densities in a cavity.
#
# cav: RCLED object
#
# source: function giving free space source amplitudes for the current kt.
#
# The resulting densities can be calculated per different quantities:
#
#   density='per_solid_angle': returns the power density per unit solid angle: 
#     Total power = Int(P(theta) 2 pi sin(theta) d theta)
#
#   density='per_kt': returns the power density per unit kt: 
#     Total power = Int(P(kt) d kt)
#
#   density='per_kx_ky': returns the power density per unit kx ky: 
#     Total power = Int(P(kx_ky) d kx d ky)
#
# The function returns P_source_top, P_source_bot, P_top (power emitted into
# top exit medium), P_sub (power emitted into substrate), P_out (power
# emitted into bottom exit medium under the substrate).
# E.g. the total outcoupling efficiency for this kt is
# P_out / (P_source_top + P_source_bot)
#
# Note that this does not account for power coming from evanescent waves
# in the source layer. A way around is the good old high index reference
# layer from Hans Deneve. This scales the power fluxes, but the extraction
# efficiency is unmodified. However, if you want to calculate the
# extraction efficiency for a an average orientation, you'll need to
# tinker with certain correction factors, see the function 'summary'.
#
#############################################################################

def P(cav, source, kt, density='per_solid_angle'):

  # Auxiliary variables.

  top = cav.top
  bot = cav.bot
  sub = cav.sub

  n_top_inc = top.inc().core().n()
  n_top_ext = top.ext().core().n()
  
  n_bot_inc = bot.inc().core().n()
  n_bot_ext = bot.ext().core().n()

  n_sub_inc = sub.inc().core().n()
  n_sub_ext = sub.ext().core().n()
    
  k0 = 2.*pi/get_lambda()

  kz_top_inc_2 = (k0*n_top_inc)**2 - kt**2
  kz_top_ext_2 = (k0*n_top_ext)**2 - kt**2
  
  kz_bot_inc_2 = (k0*n_bot_inc)**2 - kt**2
  kz_bot_ext_2 = (k0*n_bot_ext)**2 - kt**2  

  kz_sub_inc_2 = (k0*n_sub_inc)**2 - kt**2
  kz_sub_ext_2 = (k0*n_sub_ext)**2 - kt**2
    
  # Cavity modified radiation profile.

  bot.inc().set_kt(kt)
    
  top.calc()
  bot.calc()
  sub.calc()

  R_top = top.R12(0,0)
  R_bot = bot.R12(0,0)

  A_top_0, A_bot_0 = source(kt, cav)

  A_top = (A_top_0 + R_bot * A_bot_0) / (1 - R_bot * R_top)
  A_bot = (A_bot_0 + R_top * A_top_0) / (1 - R_bot * R_top)

  # Power transmission in normal direction

  T_top = abs(top.T12(0,0))**2 * top.ext().mode(0).field(Coord(0,0,0)).Sz() \
                               / top.inc().mode(0).field(Coord(0,0,0)).Sz()
  
  T_bot = abs(bot.T12(0,0))**2 * bot.ext().mode(0).field(Coord(0,0,0)).Sz() \
                               / bot.inc().mode(0).field(Coord(0,0,0)).Sz()
  
  T_sub = abs(sub.T12(0,0))**2 * sub.ext().mode(0).field(Coord(0,0,0)).Sz() \
                               / sub.inc().mode(0).field(Coord(0,0,0)).Sz()
  
  # Solid angle transformation.

  top_angle_transform = (n_top_ext * safe_sqrt(kz_top_ext_2)) / \
                        (n_top_inc * safe_sqrt(kz_top_inc_2))

  bot_angle_transform = (n_bot_ext * safe_sqrt(kz_bot_ext_2)) / \
                        (n_bot_inc * safe_sqrt(kz_bot_inc_2))
  
  sub_angle_transform = (n_sub_ext * safe_sqrt(kz_sub_ext_2)) / \
                        (n_sub_inc * safe_sqrt(kz_sub_inc_2))
    
  # Power densities per unit solid angle.
  
  P_source_top = (1.0 + 0.0j) * abs(A_top)**2 * (1 - abs(R_top)**2)
  P_source_bot = (1.0 + 0.0j) * abs(A_bot)**2 * (1 - abs(R_bot)**2)
  
  P_top = (1.0 + 0.0j) * abs(A_top)**2 * T_top * top_angle_transform
  P_sub = (1.0 + 0.0j) * abs(A_bot)**2 * T_bot * bot_angle_transform

  # Neglect Fresnel coefficient.
  # if abs(T_sub.real) > 0.01:
  #   T_sub = 1.0

  P_out = P_sub * T_sub * sub_angle_transform 
  
  # Power density per unit kt.

  if density == 'per_kt':
    
    P_source_top *= 2.*pi * kt / k0 / n_top_inc / safe_sqrt(kz_top_inc_2)
    P_source_bot *= 2.*pi * kt / k0 / n_bot_inc / safe_sqrt(kz_bot_inc_2)
    
    P_top        *= 2.*pi * kt / k0 / n_top_ext / safe_sqrt(kz_top_ext_2)
    P_sub        *= 2.*pi * kt / k0 / n_bot_ext / safe_sqrt(kz_bot_ext_2)
    
    P_out        *= 2.*pi * kt / k0 / n_sub_ext / safe_sqrt(kz_sub_ext_2)

  # Power density per unit kx ky.

  elif density == 'per_kx_ky':
    
    P_source_top *= 1.0 / k0 / n_top_inc / safe_sqrt(kz_top_inc_2)
    P_source_bot *= 1.0 / k0 / n_bot_inc / safe_sqrt(kz_bot_inc_2)
    
    P_top        *= 1.0 / k0 / n_top_ext / safe_sqrt(kz_top_ext_2)
    P_sub        *= 1.0 / k0 / n_bot_ext / safe_sqrt(kz_bot_ext_2)
  
    P_out        *= 1.0 / k0 / n_sub_ext / safe_sqrt(kz_sub_ext_2)

  elif density != 'per_solid_angle':
    print("Unknown density. Assuming per solid angle.")
  
  return P_source_top.real, P_source_bot.real, \
         abs(P_top.real), abs(P_sub.real), abs(P_out.real)



#############################################################################
#
# Integrate a peaked function f(x,args) with peak locations given by
# 'points'.
#
# In the interval points[i]-eps, points[i]+eps, a Romberg routine with
# 2**N intervals is used.
#
# In the other intervals we use quadpack's adaptive integration routines:
#
#    http://www.scipy.org/documentation/apidocs/scipy 
#                        /scipy.integrate.quadpack.html
#
# Note: as long as there is some loss in the cavity, the peaks are handled
# fine by just the adaptive integral, i.e. without providing a 'points'
# vector.
#
# TODO: simplify this, as quad itself also has an option to pass problem
# points to the function.
#
#############################################################################

def integrate_peaked_function(f, x0, x1, args=None, points=None,eps=0.1,N=10):

  import scipy.integrate.quadpack
  import scipy.integrate
  quad = scipy.integrate.quadpack.quad
  romb = scipy.integrate.romb
  
  limit = 250
  
  if points == None or points == []:
    return quad(f, x0, x1, args, limit=limit)[0]

  points.sort()

  if points[0] < x0 or points[-1] > x1:
    raise Exception("Bad peak values in integrate_peaked_function")

  result = 0
  
  # Until first peak.

  result += quad(f, x0, points[0]-eps, args, limit=limit)[0]

  # Rest of the intervals.

  for i in range(len(points)):

    # Peak.

    x0_i = points[i] - eps
    x1_i = points[i] + eps
  
    y = []
    dx = (x1_i - x0_i) / 2**N
    for j in range(2**N+1):
      y.append(f(x0_i + j*dx, args))

    result += romb(y, dx)

    # Well behaved part.

    x0_i = points[i] + eps

    if i == len(points)-1:
      x1_i = x1
    else:
      x1_i = points[i+1]-eps      
    
    result += quad(f, x0_i, x1_i, args, limit=limit)[0]

  return result



#############################################################################
#
# Calculate power fluxes.
#
#############################################################################

def calc_power_fluxes(cav, source, waveguide=None):

  # Shortcut for TE vertical.

  if source == vertical and get_polarisation() == TE:
    return 0.0, 0.0, 0.0, 0.0, 0.0

  # Other cases.

  k0 = 2.*pi/get_lambda().real
  
  kt_end_in  = k0 * cav.top.inc().core().n().real-.0001
  kt_end_top = k0 * cav.top.ext().core().n().real-.0001
  kt_end_bot = k0 * cav.bot.ext().core().n().real-.0001
  kt_end_sub = k0 * cav.sub.ext().core().n().real-.0001

  def f(kt, index):
    result = P(cav, source, kt, density='per_kt')
    return result[index]

  points = []

  if waveguide:
      
    waveguide.calc()
  
    for i in range(waveguide.N()):
      if abs(waveguide.mode(i).kz().real) > abs(waveguide.mode(i).kz().imag):
        points.append(waveguide.mode(i).kz().real)
        
  P_source_top = integrate_peaked_function(f, 0, kt_end_in,  0, points)
  P_source_bot = integrate_peaked_function(f, 0, kt_end_in,  1, points)
  P_top        = integrate_peaked_function(f, 0, kt_end_top, 2)
  P_sub        = integrate_peaked_function(f, 0, kt_end_bot, 3)
  P_out        = integrate_peaked_function(f, 0, kt_end_bot, 4)
  
  return P_source_top, P_source_bot, P_top, P_sub, P_out



#############################################################################
#
# Object holding all the results for different source orientations.
#
#############################################################################

class Result:

  def __init__(self, results):

    self.results = results

    P_source, P_sub, P_out, norm = 0.0, 0.0, 0.0, 0.0

    for source in results.keys():
    
      res = self.results[source]
    
      P_source += res.weight * res.P_source
      P_sub    += res.weight * res.P_sub
      P_out    += res.weight * res.P_out

      norm += res.weight

    self.P_source = P_source/norm
    self.eta_sub  = P_sub/P_source
    self.eta_out  = P_out/P_source

  def __getitem__(self, item):
    return self.results.__getitem__(item)

  def report(self):

    for source in self.results.keys():
      
      res = self.results[source]

      gen     = res.P_source
      eta_sub = res.eta_sub
      eta_out = res.eta_out
  
      print("")
      print("Source orientation:", source.__name__, \
            '(weight='+str(res.weight)+')'   )
      print("Emitted power                                 :",'%.3f'% gen  )
      print("Extraction efficiency to substrate            :",'%.3f'% eta_sub   )
      print("Extraction efficiency to bottom outside world :",'%.3f'% eta_out  )

    gen     = self.P_source
    eta_sub = self.eta_sub
    eta_out = self.eta_out
      
    print("")
    print("Source orientation: average" )
    print("Emitted power                                 :",'%.3f' % gen  )
    print("Extraction efficiency to substrate            :",'%.3f' % eta_sub   )
    print("Extraction efficiency to bottom outside world :",'%.3f' % eta_out  )



#############################################################################
#
# Object holding the result for a single source orientation.
#
#############################################################################

class SourceResult:

  def __init__(self, weight, P_source, P_sub, P_out):

    self.weight   = weight
    
    self.P_source = P_source
    self.P_sub    = P_sub
    self.P_out    = P_out
    
    self.eta_sub  = P_sub/P_source
    self.eta_out  = P_out/P_source



#############################################################################
#
# calc
#
#  Calculate the power flows for a set of sources which are averaged with
#  a set of weight coefficients.
#
#############################################################################

def calc(cav, sources=None, weights=None):

  # Check input. This allows you to omit specifying a weight vector,
  # rather then using [1,2] in the general case.

  if sources == None:
    sources = [vertical,horizontal]
    weights = [1,2]

  try:
    sources[0]
  except:
    sources = [sources]
    
  if weights == None:
    weights = len(sources)*[1]    
  
  if len(sources) != len(weights):
    raise Exception("Number of weights does not match number of sources.")

  # Do calculation.

  results = {}

  for source, weight in zip(sources, weights):

    # Correction factors. TODO: is this really valid for the complex case?

    if source == horizontal:
      C = cav.cor.real
    else:
      C = cav.cor.real**5   
  
    # TE.

    set_polarisation(TE)

    P_source_top_TE, P_source_bot_TE, P_top_TE, P_sub_TE, P_out_TE \
       = C*array(calc_power_fluxes(cav, source))

    # TM.

    set_polarisation(TM)

    P_source_top_TM, P_source_bot_TM, P_top_TM, P_sub_TM, P_out_TM \
       = C*array(calc_power_fluxes(cav, source))

    if source in results.keys():
      raise Exception("Error: calculated same source twice.")

    results[source] = SourceResult(weight,
                                   P_source_top_TE + P_source_bot_TE + \
                                   P_source_top_TM + P_source_bot_TM,
                                   P_sub_TE + P_sub_TM,
                                   P_out_TE + P_out_TM)

  # Print out and return results.

  r = Result(results)

  r.report()

  return r



#############################################################################
#
# Calculate angular distribution of power fluxes.
#
#############################################################################

def radiation_profile(cav, source, steps=250,
                      location='out', density='per_solid_angle',
                      show_=True, filename=None):

  if location == 'source_top' or location == 'source_bot':
    n = cav.top.inc().core().n()
  elif location == 'top':
    n = cav.top.ext().core().n() 
  elif location == 'sub':
    n = cav.bot.ext().core().n()
  elif location == 'out':
    n = cav.sub.ext().core().n() 
  else:
      raise Exception("Unknown location for radiation profile.")

  P_TE, P_TM = [], []

  theta_range = linspace(-pi/2+1e-7, pi/2-1e-7, num=steps)
  
  for theta in theta_range:

    kt = 2.*pi/get_lambda() * n * sin(theta)

    set_polarisation(TE)
    
    P_source_top, P_source_bot, P_top, P_sub, P_out = \
         P(cav, source, kt, density)
    
    if location == 'source_top':
      P_TE.append(P_source_top)
    elif location == 'source_bot':
      P_TE.append(P_source_bot)
    elif location == 'top':
      P_TE.append(P_top)
    elif location == 'sub':
      P_TE.append(P_sub)
    elif location == 'out':
      P_TE.append(P_out)

    set_polarisation(TM)
    
    P_source_top, P_source_bot, P_top, P_sub, P_out = \
         P(cav, source, kt, density)
    
    if location == 'source_top':
      P_TM.append(P_source_top)
    elif location == 'source_bot':
      P_TM.append(P_source_bot)
    elif location == 'top':
      P_TM.append(P_top)
    elif location == 'sub':
      P_TM.append(P_sub)
    elif location == 'out':
      P_TM.append(P_out)

  # Make plot.

  if show_ == True or filename != None:
    
    figure(1)
    polar(theta_range, P_TE)
    title("TE radiation density")
    if filename != None:
      savefig(filename+'_TE.png')

    figure(2)
    polar(theta_range, P_TM)
    title("TM radiation density")
    if filename != None:
      savefig(filename+'_TM.png')
      
    if show_ == True:
      show()
    
  # Return results.
      
  return theta_range, P_TE, P_TM



#############################################################################
#
# Uniform
#
#  Convenience function which creates Planars depending on
#  how it is called:
#
#  Uniform(air), Uniform(air,.5), Uniform(1), Uniform(1, .5)
#
#############################################################################

cache = []

def Uniform(material, thickness=None):
  
  global cache

  if not isinstance(material, Material):
    material = Material(material)
    cache.append(material)

  planar = Planar(material)

  cache.append(planar)

  if thickness != None:
    return planar(thickness)
  else:
    return planar



#############################################################################
#
# Main
#
#############################################################################

if __name__ == "__main__":

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


  
