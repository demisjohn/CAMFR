#! /usr/bin/env python

##############################################################################
#
# Planar TE test
#
##############################################################################

from camfr import *
from cmath import *

import unittest, eps

class planarTE(unittest.TestCase):
    def testRT(self):

        """Planar TE"""

        set_lambda(1)
        set_N(1)

        set_polarisation(TE) # needs to be set before Planars are created.

        GaAs_m = Material(3.5)
        air_m  = Material(1.0)

        GaAs = Planar(GaAs_m)
        air  = Planar( air_m)

        s = Stack(GaAs(0) + 2 * (GaAs(.1) + air(.1)) + air(0))

        GaAs.set_theta(45*pi/180.)

        s.calc()

        R = s.R12(0,0)
        T = s.T12(0,0)

        R_OK = 0.976500927535-0.215513198026j
        T_OK = -0.951371006548+0.103735346275j

        print R, "expected", R_OK
        print T, "expected", T_OK

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        self.failUnless(R_pass and T_pass)

suite = unittest.makeSuite(planarTE, 'test')        

if __name__ == "__main__":
    unittest.main()
    
