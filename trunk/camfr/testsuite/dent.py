#! /usr/bin/env python

###################################################################
#
# Dent in waveguide.
#
###################################################################

from camfr_work import *

import unittest, eps

class dent(unittest.TestCase):
    def testdent(self):

        """Dent"""

        print
        print "Running dent..."
	
	set_lambda(1.55)
	set_N(40)
        set_mode_surplus(2)

        core = Material(3.45)
        clad = Material(2)
        air  = Material(1)

	spot_size = 0.001

        set_upper_PML(-0.25)
        set_lower_PML(-0.25)

        wg = Slab(clad(6) + core(0.22) + clad(6))
        gap = Slab(air(wg.width()))
        notch = Slab(clad(6) + air(spot_size) \
                      + core(0.22-spot_size) + clad(6))

        s = Stack(wg(0) + notch(spot_size) + gap(0.28) + wg(0))

        s.calc()

        R = s.R12(0,0)
	R_OK = 0.675552207942+0.392863210495j
	
        print R, "expected", R_OK
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        set_upper_PML(0)
        set_lower_PML(0)
           
        self.failUnless(R_pass)

suite = unittest.makeSuite(dent, 'test')

if __name__ == "__main__":
    unittest.main()
