#! /usr/bin/env python

#############################################################################
#
# Polaritonic material.
#
#############################################################################

from camfr import *
from cmath import *

import unittest, eps

class polariton(unittest.TestCase):
    def testpolariton(self):

        """Polariton"""

        print('')
        print("Running polariton...")

        set_N(40)
        set_polarisation(TM)
        set_mode_surplus(2)

        # Geometry parameters.

        a = .600        # period
        r = a*.250/2.0  # rods radius

        f = 0.6
        set_lambda(a/f)

        # Materials.

        air = Material(1.0)

        epsinf = 5.1
        fL = 1.0
        fT = 0.4

        eps_r = epsinf*(f*f-fL*fL)/(f*f-fT*fT)
        polmat = Material(sqrt(eps_r))

        # Calculate band diagram

        slab1 = Slab(air(0.5*a - r) + polmat(r))
        slab2 = Slab(air(0.5*a))

        s = BlochStack(slab1(2*r) + slab2(a-2*r))
        s.calc()

        kz = s.mode(0).kz()
        kz_OK = 3.80047194128-3.70167318519e-05j
        
        print(kz, "expected", kz_OK)
        kz_pass = abs((kz - kz_OK) / kz_OK) < eps.testing_eps

        free_tmps()
        set_polarisation(TE)
           
        self.failUnless(kz_pass)

suite = unittest.makeSuite(polariton, 'test')

if __name__ == "__main__":
    unittest.main()
