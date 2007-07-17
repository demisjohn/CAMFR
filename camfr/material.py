#! /usr/bin/env python

from camfr import *

##############################################################################
#
# A simple factory class which creates a dispersive material based on a
# three-column text file in the format:
#
#   lambda (in nm)  n_real n_imag (pos means loss)
#
# When calling this object, a new Material is created with the correct
# refractive index based on the current values of lambda in micron.
#
##############################################################################

class Dispersive_Material_Factory:

  def __init__(self, filename):

    import scipy.interpolate.interpolate
    interp1d = scipy.interpolate.interpolate.interp1d

    f = file(filename)

    wavelength = []
    n = []

    for line in f:
      a,b,c  = line.split()
      wavelength.append(float(a)/1000.) # Convert to micron.
      n.append(float(b) - float(c)*1j )

    self.interpolate = interp1d(wavelength, n)

  def __call__(self):

    return Material(self.interpolate(get_lambda().real)[0])



#############################################################################
#
# Illustration of how to write a factory class for a material with
# dispersion given by a formula, in this case ZnS, from DeVore, JOSA 41,
# 416 (1951)
#
#############################################################################

class ZnS_Factory:

  def __call__(self):

    eps = 5.164 + 0.1208 / (get_lambda().real **2 - 0.27055**2)
    return Material(sqrt(eps))



#############################################################################
#
# Example
#
#############################################################################

if __name__ == "__main__":

  SiO2_factory = Dispersive_Material_Factory('SiO2.txt')

  for wavelength in arange(1,2,0.1):
    set_lambda(wavelength)
    m = SiO2_factory()
    print m.n()



