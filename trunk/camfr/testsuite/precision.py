#!/usr/bin/env python

###################################################################
#
# Test for high precision.
#
###################################################################

from camfr_work import *

import unittest, eps

class precision(unittest.TestCase):
    def testprecision(self):

        """Precision"""

        print
        print "Running precision..."
        
        set_N(20)
        set_mode_surplus(3.0)
        set_polarisation(TE)
        set_lambda(.98)

        set_precision(1000)

        m = Material(3.0)
        a = Material(1.0)

        waveguide  = Slab(a(3.0) + m(0.5) + a(3.0))
        space = Slab(a(6.5))

        s = Stack(waveguide(1) + space(1))
        s.calc()

        R = s.R12(0,0)
        R_OK = 0.449374475028+0.415450001922j

        print R, "expected", R_OK
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()
        set_precision(100)
 
        self.failUnless(R_pass)

suite = unittest.makeSuite(precision, 'test')

if __name__ == "__main__":
    unittest.main()
          
