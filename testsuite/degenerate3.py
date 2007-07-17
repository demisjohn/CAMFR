#! /usr/bin/env python

###################################################################
#
# Test for uncoupled waveguides.
#
###################################################################

from camfr import *

import unittest, eps

class degenerate3(unittest.TestCase):
    def testdegenerate3(self):

        """Degenerate 3"""

        print
        print "Running degenerate 3..."
        
        set_lambda(1.5)
        set_N(2)

        core=Material(3.2)
        clad=Material(1.0)

        set_unstable_exp_threshold(1e-6)

        s1 = Slab(clad(1)+core(0.3)+clad(5.764))
        s2 = Slab(clad(1)+core(0.3)+clad(4.464)+core(0.3)+clad(1))

        s = Stack(s1(0)+s2(0))

        s.calc()

        R = s.R12(0,0)
        R_OK = 0.0

        print R, "expected", R_OK
        R_pass = abs(R - R_OK) < eps.testing_eps

        free_tmps()

        set_unstable_exp_threshold(1e-12)
        
        self.failUnless(R_pass)

suite = unittest.makeSuite(degenerate3, 'test')

if __name__ == "__main__":
    unittest.main()
          
