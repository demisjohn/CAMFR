#! /usr/bin/env python
  
from camfr import *
from numpy.linalg import *
from cmath import *
import sys

warned = False

# Temporary hack to allow looping over kx0 and ky0.

set_always_recalculate(True)



#############################################################################
#
# V2.0, 20070706
#
# Library to simulate RCLEDs which have gratings in the cavity.
#
# Possibilities for improvement:
#
#  write a quad routine which can operate on vectors
#
#  more freedom of types of symmetry to use
#
#  Fix Z multiplication for TE/TM at the C++ level rather than at the
#  Python level
#
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
# for the total power then. See summary() function.
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
    
    # Create a Stack composed of the flipped bottom stack followed by the
    # top stack.
    
    self.bot_flipped     = FlippedScatterer(self.bot.scatterer())
    self.bot_flipped_top = Stack(self.bot_flipped + top)
    
    # Consistency check.

    if self.top.inc().core().n() != self.bot.inc().core().n():
      raise Exception("Incidence media of top and bottom stack don't match.")

    if self.bot.ext().core().n() != self.sub.inc().core().n():
      raise Exception("Exit medium of bottom doesn't match inc medium of sub.")

  def calc(self, sources=None, weights=None, steps=30, symmetric=False,
           single_pass_substrate=True):
    return calc(self, sources, weights, steps, symmetric, True,
                single_pass_substrate)

  def radiation_profile(self, source, steps=30,
                        location='out', density='per_solid_angle',
                        show=True, fname=None):
    return radiation_profile(self,source,steps,location,density,show,fname)



  ###########################################################################
  #
  # Gives the minimal orders needed to span all k-space.
  #
  ###########################################################################

  def minimum_orders(self):
    
    if get_lambda() == 0:
      raise Exception("Lambda needs to be set.")

    Kx = 2.*pi / self.top.inc().width()
    Ky = 2.*pi / self.top.inc().height()

    k0 = 2.*pi/get_lambda().real

    max_k_needed = k0*self.top.inc().core().n().real

    return int(abs(max_k_needed - Kx/2.0)/Kx + 1), \
           int(abs(max_k_needed - Ky/2.0)/Ky + 1)



#############################################################################
#
# Source terms for a dipole with a given orientation, which can be either
#
#   'vertical'     : perpendicular to layer interfaces
#   'horizontal_x' : parallel to interface, along x
#   'horizontal_y' : parallel to interface, along y
#
# The amplitude represent power densities per unit dx dy.
#
# When source() is called, it creates two source vectors containing the
# containing the upward and downward amplitudes for all directions coupled
# by the Bragg condition.
#
#############################################################################

class Source:

  def __init__(self, cav, orientation):

    self.wg = cav.top.inc()
    self.orientation = orientation



  def _amplitude(self, kx, ky, polarisation):

    global warned
    
    n   = self.wg.core().n()
    k0  = 2.*pi/get_lambda()
    kz2 = (k0*n)**2 - kx**2 - ky**2
    kz  = safe_sqrt(kz2)
    kt  = safe_sqrt(kx**2 + ky**2)
    phi = atan2(ky.real, kx.real)

    Z = sqrt(self.wg.core().mu() / self.wg.core().eps())
    
    if abs(kz2) < 1e-12:
      if warned == False:
        print kx, ky, "Warning: kz=0: close to cut-off."
        warned = True
      kz = 1e-12

    #
    # Vertical dipole.
    #
    
    if self.orientation == "vertical":
      
      if polarisation == TM:
        
        A  = sqrt(3./8./pi) * kt/k0/n # per unit solid angle
        A /= k0 * n * kz              # per unit dkx dky
        return A, -A

      else:
        
        return 0.0, 0.0
    
    #
    # Horizontal dipole along x.
    #
    
    if self.orientation == "horizontal_x":
       
      if polarisation == TE:

        A = sqrt(3./8./pi) * cos(phi) # per unit solid angle
        A /= k0 * n * kz              # per unit dkx dky
        return A*Z, A*Z

      else:

        A  = sqrt(3./8./pi) * kz/k0/n * sin(phi) # per unit solid angle
        A /= k0 * n * kz                         # per unit dkx dky
        return A, A
    
    #
    # Horizontal dipole along y.
    #

    if self.orientation == "horizontal_y":
      
      if polarisation == TE:

        A = -sqrt(3./8./pi) * sin(phi) # per unit solid angle
        A /= k0 * n * kz               # per unit dkx dky
        return A*Z, A*Z

      else:

        A  = sqrt(3./8./pi) * kz/k0/n * cos(phi) # per unit solid angle
        A /= k0 * n * kz                         # per unit dkx dky
        return A, A

    else:
      raise Exception("Unknown orientation.")



  def __call__(self, kx0, ky0):
    
    self.wg.set_kx0_ky0(kx0, ky0)
    self.wg.calc()

    A_0_top = array(zeros(self.wg.N()), complex)
    A_0_bot = array(zeros(self.wg.N()), complex)
    
    for i in range(self.wg.N()):
      
      m = self.wg.mode(i)
      
      A_0_top[i], A_0_bot[i] = self._amplitude(m.get_kx(),m.get_ky(),m.pol())
    
    return A_0_top, A_0_bot



#############################################################################
#
# Convenience functions to create the same interface as in the planar case.
#
#############################################################################

def vertical(cav, kx0, ky0):
  return Source(cav, "vertical")(kx0, ky0)

def horizontal_x(cav, kx0, ky0):
  return Source(cav, "horizontal_x")(kx0, ky0)

def horizontal_y(cav, kx0, ky0):
  return Source(cav, "horizontal_y")(kx0, ky0)



#############################################################################
#
# Safe matrix inversion, uses svd if matrix is singular.
#
#############################################################################

def inv_safe(A):

  try:

    return inv(A)

  except:

    print "Singularities dectected in matrix inversion."

    # A = U * sigma * Vh

    U, sigma, Vh = linalg.svd(A)

    # inv_A = V * 1/sigma * Uh, but set 1/0 in sigma to 0.

    inv_sigma = where(abs(sigma) <= 1e-11, 0, 1/sigma) 
    
    return dot(Vh.conj().T, dot(diag(inv_sigma), U.conj().T))



#############################################################################
#
# Calculate extracted field through a cavity, but only for waves
# coupled to kx0 and ky0 through the Bragg condition.
#
# cav: RCLED cavity
#
# sources: single source field distribution, or a list of sources
#
# if add_all == True, then all contributions for all kx and ky for this kx0
# and ky0 will be summed.
#
#############################################################################

def P(cav, sources, kx0, ky0, single_pass_substrate=True, add_all=True,
      density='per_kx_ky'):

  # If we only got passed a single source, turn it into a list.

  try:
    sources[0]
  except:
    sources = [sources]

  # Calculate scattering matrices.
  
  cav.top.inc().set_kx0_ky0(kx0, ky0)
  
  cav.top.calc()
  cav.bot.calc()

  R_top = cav.top.R12()
  R_bot = cav.bot.R12()

  N = cav.top.inc().N()

  inv_1 = inv_safe(identity(N) - dot(R_bot, R_top)) 
  inv_2 = inv_safe(identity(N) - dot(R_top, R_bot))

  # Loop over the different sources.

  results = []
  
  for source in sources:
  
    # Cavity modified radiation profile.

    A_top_0, A_bot_0 = source(cav, kx0, ky0)

    A_top_up   = dot(inv_1,A_top_0 + dot(R_bot, A_bot_0))
    A_top_down = dot(R_top,A_top_up)

    A_bot_down = dot(inv_2,A_bot_0 + dot(R_top, A_top_0))
    A_bot_up   = dot(R_bot,A_bot_down)

    # All C++ implementation.

    # TODO: investigate why this all C++ technique takes twice as long
    # Also see if always_recalculate(True) can be gotten rid of, allowing
    # to cache the matrices which are independent from the source also in C++.

    #c = Cavity(cav.bot, cav.top)
    #c.set_source(A_top_0, A_bot_0)
    
    #A_top_up   = cav.top.inc_field()
    #A_top_down = cav.top.refl_field()
    
    #A_bot_down = cav.bot.inc_field()
    #A_bot_up   = cav.bot.refl_field()
     
    # Net power radiated by the source (take fw and bw into account).

    P_source_top = abs(A_top_up  )**2 - abs(A_top_down)**2
    P_source_bot = abs(A_bot_down)**2 - abs(A_bot_up  )**2

    for i in range(N):
      p = cav.top.inc().mode(i).field(Coord(0,0,0)).Sz().real
      P_source_top[i] *= p
      P_source_bot[i] *= p
    
    #for i in range(N): # Todo: hoist out of loop
    #  m = cav.top.inc().mode(i)
    #  P_source_top[i] = (abs(A_top_up  [i])**2 - abs(A_top_down[i])**2) \
    #                       * m.field(Coord(0,0,0)).Sz().real
    #  P_source_bot[i] = (abs(A_bot_down[i])**2 - abs(A_bot_up  [i])**2) \
    #                       * m.field(Coord(0,0,0)).Sz().real
    
    # Power transmitted into substrate.
    
    P_sub = abs(dot(cav.bot.T12(), A_bot_down))**2

    for i in range(N):
      m = cav.bot.ext().mode(i)
      P_sub[i] *= m.field(Coord(0,0,0)).Sz().real

    # Power transmitted to bottom outside world.

    if single_pass_substrate:

      cav.sub.calc()

      P_out = dot(cav.sub.T12_power(), P_sub).real

    else:
              
      cav.bot_flipped_top.calc()
      cav.sub.calc()
 
      R_cav = cav.bot_flipped_top.R12_power().real
      R_sub = cav.sub.R12_power().real
      T_sub = cav.sub.T12_power().real

      # Remove evanescent modes from the picture.
        
      propagating = []
        
      for i in range(N):
        if abs(cav.sub.inc().mode(i).field(Coord(0,0,0)).Sz().real) > 1e-8:
          propagating.append(i)
        
      P_sub_  =    take(P_sub,  propagating)
      R_cav = take(take(R_cav,propagating,0),propagating, 1)
      R_sub = take(take(R_sub,propagating,0),propagating, 1)
      T_sub = take(take(T_sub,propagating,0),propagating, 1)

      # Calculate transmitted power density.
        
      if len(propagating) > 0: # LAPACK crashes on empty matrix.
        inv_  = inv_safe(identity(len(propagating)) - dot(R_cav, R_sub))
        P_out = dot(T_sub, dot(inv_, P_sub_))
      else:
        P_out = [0.0]

      # Note that P_out can now have a different size than P_source. When
      # just calculating total power, this is no problem. It is however an
      # issue for field profiles.

    # Transform to solid angle. Only used for field profiles. Not valid for
    # multipass substrate.
    
    if density == 'per_solid_angle':

      k0 = 2.*pi/get_lambda()

      n_top_inc = cav.top.inc().core().n()
      n_bot_ext = cav.bot.ext().core().n()  
      n_sub_ext = cav.sub.ext().core().n()

      for i in range(N):
        c = cav.top.inc().mode(i).kz()
        P_source_top[i] *= k0 * n_top_inc * c
        P_source_bot[i] *= k0 * n_top_inc * c
        
      for i in range(N):
        c = cav.bot.ext().mode(i).kz()    
        P_sub[i] *= k0 * n_bot_ext * c
        
      for i in range(N):
        c = cav.sub.ext().mode(i).kz()
        P_out[i] *= k0 * n_sub_ext * c
              
    # Correction factor taking the presence of the reference layer
    # into account. Gives correct power balance between vertical and
    # horizontal dipoles.
    # TODO: is this really valid for the lossy case?

    if source == vertical:
      C = (cav.cor**5).real 
    else:
      C = (cav.cor).real
      
    if add_all == True: # TODO: refactor
      results += [C*sum(P_source_top),
                  C*sum(P_source_bot),
                  0.0, # Top emission not calculated to save time.
                  C*sum(P_sub), C*sum(P_out)]
    else:
      if single_pass_substrate == False:
        print "Warning: fields matrices not yet implemented",
        print "for multipass substrate. P_out can be inconsistent."
      results += [C*P_source_top, C*P_source_bot, None, C*P_sub, C*P_out]
  
  return results



#############################################################################
#
# Calculate radiation profile.
#
#############################################################################

def radiation_profile(cav, source, steps=30,
                      location='out', density='per_solid_angle',
                      show_=True, fname=None):

  if location == 'source_top':
    n_mat = cav.top.inc().core().n()
    index = 0
  elif location == 'source_bot':
    n_mat = cav.bot.inc().core().n()
    index = 1
  elif location == 'top':
    raise Exception("Top radiation not calculated to save time.")
  elif location == 'sub':
    n_mat = cav.bot.ext().core().n()
    index = 3
  elif location == 'out':
    n_mat = cav.sub.ext().core().n()
    index = 4   
  else:
    raise Exception("Unknown location for radiation profile.")

  # Integration limits. In case of a symmetric grating, these could be
  # reduced.

  Kx = 2.*pi / cav.top.inc().width()
  Ky = 2.*pi / cav.top.inc().height()

  kx0, kx1 = -Kx/2., Kx/2.
  ky0, ky1 = -Ky/2., Ky/2.
  
  k0 = 2.*pi/get_lambda().real

  Mx = get_fourier_orders_x()
  My = get_fourier_orders_y()
  
  max_k_spanned = Kx/2. + Mx*Kx
  max_k_needed  = k0*cav.top.inc().core().n().real

  if max_k_spanned < max_k_needed:
    Mx_required, My_required = cav.minimum_orders()
    
    print "Warning: not all of k-space is covered."
    print "Suggested orders: Mx:", Mx_required, ", My:", My_required 

  # Set up the data structures.
  
  kx0_range = linspace(kx0, kx1, num=steps)
  ky0_range = linspace(ky0, ky1, num=steps)
  
  # Create list of all kx and ky for a given kx0 and ky0.

  kx_list = []
  for kx0 in kx0_range:
    for n in range(-Mx, Mx+1, 1):
      kx = kx0 + n*Kx
      if len(where(abs(array(kx_list) - kx) < 1e-12)[0]) == 0: # Not present.
        kx_list.append(kx)
  kx_list.sort()
  
  ky_list = []
  for ky0 in ky0_range:
    for m in range(-My, My+1, 1):
      ky = ky0 + m*Ky
      if len(where(abs(array(ky_list) - ky) < 1e-12)[0]) == 0: # Not present.
        ky_list.append(ky)
  ky_list.sort()

  M = len(kx_list)
  
  P_TE = zeros((M,M), float)
  P_TM = zeros((M,M), float)

  # Calculate power.

  X, Y = meshgrid(kx_list, ky_list)

  wg = cav.bot.ext()

  progress = 0.

  for kx0 in kx0_range:

    progress += 1

    #print 'Progress: %.1f'% (progress / len(kx0_range) * 100), '%'
    
    for ky0 in ky0_range:
      
      result = P(cav, source, kx0, ky0, add_all=False, density=density)[index]
      
      for i in range(len(result)):

        i_x = where(abs(array(kx_list) - wg.mode(i).get_kx().real) < 1e-12)[0]
        i_y = where(abs(array(ky_list) - wg.mode(i).get_ky().real) < 1e-12)[0]

        # Sanity checks, shouldn't happen.

        if (len(i_x) > 1) or (len(i_y) > 1):
          print "Warning: multiple indices", i_x, i_y

        #if wg.mode(i).pol() == TE and abs(P_TE[i_x,i_y]) > 1e-12:
          
          # Points at the edges of the BZ can be visited more than once, but
          # they should always give the same result.
          # If they don't, that usually disappears when taking more orders.
          
          #if abs(P_TE[i_x,i_y] - result[i]) > 1e-4:
          #  print "BZ consistency problem", Kx/2., \
          #        wg.mode(i).get_kx().real, wg.mode(i).get_ky().real, \
          #        abs(P_TE[i_x,i_y] - result[i])

        if wg.mode(i).pol() == TE:
          P_TE[i_x,i_y] = result[i]
        else:
          P_TM[i_x,i_y] = result[i]

  # Plot.

  if show_ == True or fname != None:

    X, Y = meshgrid(kx_list, ky_list)

    M = 1.4*(2.*pi/get_lambda() * n_mat).real

    figure(1)
    pcolor(X, Y, P_TE.T, shading="flat", cmap=cm.hot)

    axis('equal') # Only works when displaying the full data
    v = axis()    # For cropping, we need to save the aspect ratio.
    ratio = v[1]/v[3]

    axis([-M*ratio, M*ratio, -M, M])
    xlabel("kx")
    ylabel("ky")    
    title("TE radiation density")
    if fname != None:
      savefig(fname+'_TE.png')

    figure(2)
    pcolor(X, Y, P_TM.T, shading="flat", cmap=cm.hot)
    axis([-M*ratio, M*ratio, -M, M])
    xlabel("kx")
    ylabel("ky")  
    title("TM radiation density")
    if fname != None:
      savefig(fname+'_TM.png')
      
    if show_ == True:
      show()

  # Return results.
  
  return array(kx_list), P_TE, P_TM



#############################################################################
#
# Integrate a 2D function based on an equidistant grid of samples.
#
# Not so accurate in the case of strong peaks, but has the advantage that
# it can work on functions which return a vector.
#
#############################################################################

def integrate_2D(f, x0, x1, dx, y0, y1, dy):
  
  import scipy.integrate
  simps = scipy.integrate.simps

  xrange = arange(x0,x1+dx/2.,dx)
  yrange = arange(y0,y1+dy/2.,dy)

  integrated_over_y = []
  for ix in range(len(xrange)):
    f_i = []
    for iy in range(len(yrange)):
      f_i.append(f(xrange[ix],yrange[iy]))
    integrated_over_y.append(simps(f_i, yrange, axis=0))

  return simps(integrated_over_y, xrange, axis=0)



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
  
      print
      print "Source orientation:", source.__name__, \
            '(weight='+str(res.weight)+')'
      print "Emitted power                                 :",'%.3f'% gen
      print "Extraction efficiency to substrate            :",'%.3f'% eta_sub 
      print "Extraction efficiency to bottom outside world :",'%.3f'% eta_out

    gen     = self.P_source
    eta_sub = self.eta_sub
    eta_out = self.eta_out
      
    print
    print "Source orientation: average"
    print "Emitted power                                 :",'%.3f' % gen
    print "Extraction efficiency to substrate            :",'%.3f' % eta_sub 
    print "Extraction efficiency to bottom outside world :",'%.3f' % eta_out

    sys.stdout.flush()



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
# Calculate extraction efficiency.
# We only need to integrate over the first Brillouin zone.
#
# If steps == 0, we use the slower, but more accurate quad routines.
# Unfortunately, this routine can only handle scalars so far, so it requires
# several independent integration runs.
#
# If symmetric==True, we only integrate over a quarter of the Brillouin
# zone, which is good for a dipole along the axes in a structure with
# two mirror axes.
#
#############################################################################

def calc(cav, sources=None, weights=None, steps=30, symmetric=False,
         coherent=True, single_pass_substrate=True):
  
  # Check input. This allows you to omit specifying a weight vector,
  # rather then using [1, 1, 1] in the general case.

  if sources == None:
    sources = [vertical, horizontal_x, horizontal_y]
    weights = [1, 1, 1]

  try:
    sources[0]
  except:
    sources = [sources]
    
  if weights == None:
    weights = len(sources)*[1]    
  
  if len(sources) != len(weights):
    raise Exception("Number of weights does not match number of sources.")
  
  # Integration limits. In case of a symmetric grating, these can be
  # reduced. (TODO: clean code duplication)

  Kx = 2.*pi / cav.top.inc().width()
  Ky = 2.*pi / cav.top.inc().height()

  kx0, kx1 = -Kx/2., Kx/2.
  ky0, ky1 = -Ky/2., Ky/2.
  
  k0 = 2.*pi/get_lambda().real
  
  max_k_spanned = Kx/2. + get_fourier_orders_x()*Kx
  max_k_needed  = k0*cav.top.inc().core().n().real

  if max_k_spanned < max_k_needed:
    Mx_required, My_required = cav.minimum_orders()
    
    print "Warning: not all of k-space is covered."
    print "Suggested orders: Mx:", Mx_required, ", My:", My_required
    
  results = {}

  # If steps is set, use an equidistant grid of samples.

  if steps:
 
    # Calculation of power flows.

    def f(kx0, ky0):
      return [i for i in P(cav, sources, kx0, ky0,
                           single_pass_substrate=single_pass_substrate)]

    if symmetric:
      kx0 = ky0 = 0.0
  
    fluxes = integrate_2D(f, kx0, kx1, (kx1-kx0)/steps, \
                             ky0, ky1, (ky1-ky0)/steps)

    # Collect results.

    for i in range(len(fluxes)/5):
      gen, sub, out = fluxes[5*i]+fluxes[5*i+1],fluxes[5*i+3],fluxes[5*i+4]

      if sources[i] in results.keys():
        raise Exception("Error: calculated same source twice.")

      results[sources[i]] = SourceResult(weights[i], gen, sub, out)

    # Print out and return results.

    r = Result(results)

    r.report()

    return r

  # Else, use the slower but more accurate quad routines.
  # Note: impossibly slow, uses old interface convention.
  
  else:

    import scipy.integrate.quadpack
    dblquad = scipy.integrate.quadpack.dblquad 

    for source in sources:

      # TODO: make symmetries dependent on switch.

      if symmetric:

        gen = dblquad(lambda kx, ky : P(cav,source,kx,ky,
                             single_pass_substrate=single_pass_substrate)[0],
                      0, ky1,
                      lambda kx : 0,
                      lambda kx : sqrt((Ky/2.)**2-kx*2).real)[0]

      else:

        gen = dblquad(lambda kx, ky : P(cav,source,kx,ky,
                             single_pass_substrate=single_pass_substrate)[0],
                      ky0, ky1,
                      lambda kx : -sqrt((Ky/2.)**2-kx*2).real,
                      lambda kx :  sqrt((Ky/2.)**2-kx*2).real)[0]
         
      print
      print "Source orientation:", source.__name__
      print "Generated power                              :", '%.3f' % gen
      sys.stdout.flush()

      sub = 0.0 # Let's not calculate 'sub' to save time.
      
      if symmetric:
        out = dblquad(lambda kx, ky : P(cav,source,kx,ky,
                             single_pass_substrate=single_pass_substrate)[2],
                      0, ky1,
                      lambda ky : 0,
                      lambda ky : sqrt((Ky/2.)**2-ky*2).real)[0]
      else:
        out = dblquad(lambda kx, ky : P(cav,source,kx,ky,
                             single_pass_substrate=single_pass_substrate)[2],
                      ky0, ky1,
                      lambda ky : -sqrt((Ky/2.)**2-ky*2).real,
                      lambda ky :  sqrt((Ky/2.)**2-ky*2).real)[0]
        
      print "Outcoupled power                              :",'%.3f'% out
      print "Extraction efficiency to substrate            :",'%.3f'%(sub/gen)
      print "Extraction efficiency to bottom outside world :",'%.3f'%(out/gen)
      sys.stdout.flush()

      all_gen.append(gen)
      all_sub.append(sub)
      all_out.append(out)

  return all_gen, all_sub, all_out



#############################################################################
#
# set_period
#
#############################################################################

_Lx = None
_Ly = None

def set_period(Lx, Ly=None):
  
  global _Lx, _Ly

  if Ly == None:
    Ly = Lx
    
  _Lx = Lx 
  _Ly = Ly



#############################################################################
#
# Uniform
#
#  Convenience function which creates BlochSections or Terms depending on
#  how it is called:
#
#  Uniform(air), Uniform(air,.5), Uniform(1), Uniform(1, .5)
#
#############################################################################

cache = []

def Uniform(material, thickness=None):
  
  global cache

  if _Lx == None:
    raise Exception("period is not defined with set_period().")
  
  if not isinstance(material, Material):
    material = Material(material)
    cache.append(material)

  slab    = Slab(material(_Ly))
  section = BlochSection(slab(_Lx))  

  cache.append(slab)
  cache.append(section)

  if thickness != None:
    return section(thickness)
  else:
    return section



#############################################################################
#
# SquareGrating
#
#   Square lattice of square holes of square_m in background of surround_m.
#   The holes have side D.
#   Three different positions are available, which should be enough to
#   span the BZ in case Lx=Ly.
#
#############################################################################

def SquareGrating(D, square_m, surround_m, pos="1", thickness=None):

  global cache
  
  if _Lx == None:
    raise Exception("Period is not defined with set_period().")

  mat1 = square_m
  mat2 = surround_m

  if not isinstance(mat1, Material):
    mat1 = Material(mat1)
    cache.append(mat1)
    
  if not isinstance(mat2, Material):
    mat1 = Material(mat2)
    cache.append(mat2)

  # Define gratings. Make sure to keep mirror symmetry around origin.
  # The dipole is located in the left corner of the BlochSection.
  
  if pos == "1": # Right above hole of material mat1.
    
    # MAT1 mat2 mat2 MAT1
    # mat2 mat2 mat2 mat2
    # mat2 mat2 mat2 mat2
    # MAT1 mat2 mat2 MAT1
    
    s1 = Slab(mat1(D/2.) + mat2(_Ly-D) + mat1(D/2.))
    s2 = Slab(mat2(_Ly))  
  
    section = BlochSection(s1(D/2.) + s2(_Lx-D) + s1(D/2.))

  elif pos == "2": # Diagonally between two holes.
    
    # mat2 mat2 mat2 mat2
    # mat2 MAT1 MAT1 mat2
    # mat2 MAT1 MAT1 mat2
    # mat2 mat2 mat2 mat2
    
    s1 = Slab(mat2((_Ly-D)/2.) + mat1(D) + mat2((_Ly-D)/2.))
    s2 = Slab(mat2(_Ly))  
  
    section = BlochSection(s2((_Lx-D)/2.) + s1(D) + s2((_Lx-D)/2.))

  elif pos == "3": # Vertically between two holes.
    
    # mat2 mat2 mat2 mat2
    # MAT1 mat2 mat2 MAT1
    # MAT1 mat2 mat2 MAT1
    # mat2 mat2 mat2 mat2
    
    s1 = Slab(mat2((_Ly-D)/2.) + mat1(D) + mat2((_Ly-D)/2.))
    s2 = Slab(mat2(_Ly))  
  
    section = BlochSection(s1(D/2.) + s2(_Lx-D) + s1(D/2.))

  else:
    raise Exception("Unknown dipole position.")
  
  cache.append(s1)
  cache.append(s2)  
  cache.append(section)

  if thickness != None:
    return section(thickness)
  else:
    return section



#############################################################################
#
# SquareGratingRoundHoles
#
#   Square lattice of circular holes of core_m in background of surround_m.
#   The holes diameter D. dx and dy are the parameters used in the geometry
#   definition.
#   Three different positions are available, which should be enough to
#   span the BZ in case Lx=Ly.
#
#############################################################################

def SquareGratingRoundHoles(D, dx, dy, core_m, surround_m, pos="1",
                            thickness=None):
  
  global cache
  
  R = D/2.0
  
  if _Lx == None:
    raise Exception("Period is not defined with set_period().")

  mat1 = core_m
  mat2 = surround_m

  if not isinstance(mat1, Material):
    mat1 = Material(mat1)
    cache.append(mat1)
    
  if not isinstance(mat2, Material):
    mat1 = Material(mat2)
    cache.append(mat2)

  # Define gratings. Make sure to keep mirror symmetry around origin.
  # The dipole is located in the left corner of the BlochSection.
  
  geometry = Geometry(mat2)
  
  if pos == "1": # Right above hole of material mat1.
    
    # MAT1 mat2 mat2 MAT1
    # mat2 mat2 mat2 mat2
    # mat2 mat2 mat2 mat2
    # MAT1 mat2 mat2 MAT1

    geometry += Circle(Point(0. , 0. ), R, mat1) # lower left
    geometry += Circle(Point(0. , _Ly), R, mat1) # upper left
    geometry += Circle(Point(_Lx, 0. ), R, mat1) # lower right
    geometry += Circle(Point(_Lx, _Ly), R, mat1) # upper right

  elif pos == "2": # Diagonally between two holes.
    
    # mat2 mat2 mat2 mat2
    # mat2 MAT1 MAT1 mat2
    # mat2 MAT1 MAT1 mat2
    # mat2 mat2 mat2 mat2
    
    geometry += Circle(Point(_Lx/2.0 , _Ly/2.0), R, mat1)    
     
  elif pos == "3": # Vertically between two holes.
    
    # mat2 mat2 mat2 mat2
    # MAT1 mat2 mat2 MAT1
    # MAT1 mat2 mat2 MAT1
    # mat2 mat2 mat2 mat2
    
    geometry += Circle(Point(0.0 , _Ly/2.0), R, mat1)    
    geometry += Circle(Point(_Lx , _Ly/2.0), R, mat1)    
    
  else:
    raise Exception("Unknown dipole position.")
  
  expr = geometry.to_expression(x0=0., x1=_Lx, dx=dx,
                                y0=0., y1=_Ly, dy=dy)

  section = BlochSection(expr)
  cache.append(geometry)
  cache.append(expr)
  cache.append(section)  

  if thickness != None:
    return section(thickness)
  else:
    return section



#############################################################################
#
# Main
#
#############################################################################

if __name__ == "__main__":

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
