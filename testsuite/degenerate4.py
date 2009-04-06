#! /usr/bin/env python

###################################################################
#
# Test for uncoupled waveguides.
#
###################################################################

from camfr import *

import unittest, eps

class degenerate4(unittest.TestCase):
    def testdegenerate4(self):

        """Degenerate 4"""

        print
        print "Running degenerate 4..."


        set_lambda(1.5)
        set_N(2)

        core_width = 0.3
        clading_width = 1.5

        cladding = Material(1)
        core = Material(3.2)

        s1 = Slab(cladding(1.0) + core(core_width) +
         cladding(clading_width+ core_width + 1.0))
        s2 = Slab(cladding(1.0) + core(core_width) +
         cladding(clading_width) + core(core_width) + cladding(1.0))

        s = Stack(s1(0)+s2(0))

        s.calc()

        R = s.R12(0,0)
        R_OK = 0.0

        print R, "expected", R_OK
        R_pass = abs(R - R_OK) < eps.testing_eps

        free_tmps()
        
        self.failUnless(R_pass)

suite = unittest.makeSuite(degenerate4, 'test')

if __name__ == "__main__":
    unittest.main()
          
