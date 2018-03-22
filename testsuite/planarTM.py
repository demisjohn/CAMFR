#! /usr/bin/env python

##############################################################################
#
# Planar TM test
#
##############################################################################

from camfr import *
from cmath import *

import unittest, eps

class planarTE(unittest.TestCase):
    def testRT(self):

        """Planar TM"""

        print('')
        print("Running planar TM...")

        set_lambda(1)
        set_N(1)

        set_polarisation(TM) # needs to be set before Planars are created.

        GaAs_m = Material(3.5)
        air_m  = Material(1.0)

        GaAs = Planar(GaAs_m)
        air  = Planar( air_m)

        s = Stack(GaAs(0) + 2 * (GaAs(.1) + air(.1)) + air(0))

        GaAs.set_theta(45*pi/180.)

        s.calc()

        R = s.R12(0,0)
        T = s.T12(0,0)

        R_OK = -0.985823303384  +0.167786812667j
        T_OK =  0.00196175699271-0.000165753394309j

        print(R, "expected", R_OK)
        print(T, "expected", T_OK)

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(R_pass and T_pass)

suite = unittest.makeSuite(planarTE, 'test')        

if __name__ == "__main__":
    unittest.main()
    
