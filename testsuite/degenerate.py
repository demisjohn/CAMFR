#! /usr/bin/env python

###################################################################
#
# Test for uncoupled waveguides.
#
###################################################################

from camfr_work import *

import unittest, eps

class degenerate(unittest.TestCase):
    def testdegenerate(self):

        """Degenerate"""

        print
        print "Running degenerate..."

        set_N(10)
        set_lambda(1)
        set_polarisation(TE)

        m = Material(3.5)
        a = Material(1.0)

        s = Slab(a(2) + m(0.5) + a(30) + m(0.5) + a(2))
        s.calc()

        E0 = abs(s.mode(0).field(Coord(2.25,0,0)).E2())
        E0_OK = 19.2559587762
        
        print E0, "expected", E0_OK
        E0_pass = abs((E0 - E0_OK) / E0_OK) < eps.testing_eps
        
        E1 = abs(s.mode(1).field(Coord(s.width()-2.25,0,0)).E2())
        E1_OK = E0_OK

        print E1, "expected", E1_OK
        E1_pass = abs((E1 - E1_OK) / E1_OK) < eps.testing_eps

        free_tmps()
           
        self.failUnless(E0_pass and E1_pass)

suite = unittest.makeSuite(degenerate, 'test')

if __name__ == "__main__":
    unittest.main()
          
