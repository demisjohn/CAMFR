#! /usr/bin/env python

###################################################################
#
# Grating 3
#
###################################################################

from camfr import *

import unittest, eps

class grating3(unittest.TestCase):
    def testgrating3(self):

        """Grating 3"""

        print('')
        print("Running grating 3...")

        set_lambda(1.55)
        set_N(40)

        set_degenerate(0)

        # Define materials. 

        substrate = Material(1.444)
        guiding   = Material(3.476)
        air       = Material(1)

        # Define geometry parameters.
        
        d = 7.0

        guide_thickness = 0.220
        groove_depth = 0.050
        period_length = .550
        filling_factor = .05

        set_lower_PML(-0.4)
        set_upper_PML(-0.4)

        # Define structure.

        waveguide = Slab(substrate(d) + guiding(guide_thickness) + air(d))
        etched = Slab(substrate(d) + guiding(guide_thickness-groove_depth)
              + air(groove_depth + d))

        stack = Stack(waveguide(period_length)  \
	    + 10 * (etched(  period_length*(1-filling_factor))
                           + waveguide(period_length*filling_factor)) \
	    + waveguide(period_length))

        # Calculate.

        stack.calc()
        
        R = stack.R12(0,0)
        R_OK = 0.103410938939-0.0238576791412j
        
        print(R, "expected", R_OK)
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()
        
        set_lower_PML(0)
        set_upper_PML(0)
        set_degenerate(1)
        
        self.failUnless(R_pass)

suite = unittest.makeSuite(grating3, 'test')

if __name__ == "__main__":
    unittest.main()

