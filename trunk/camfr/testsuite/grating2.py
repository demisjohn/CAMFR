#! /usr/bin/env python

###################################################################
#
# COST268 WG2 grating.
#
###################################################################

from camfr import *

import unittest, eps

class grating2(unittest.TestCase):
    def testgrating2(self):

        """Grating 2"""

        print
        print "Running grating 2..."

        set_polarisation(TE)
        set_lambda(.84)

        d_guide = .500
        groove_depth = .625

        period_length = .430
        no_of_grooves = 20
        filling_factor = 0.5

        d_air = 4.5
        d_sub = 5

        set_N(70)
  
        # Define materials.

        substrate = Material(1.45)
        guiding   = Material(2.00)
        air       = Material(1.00)

        # Define 1D slabs.

        waveguide = Slab(substrate(d_sub) + guiding(d_guide) + air(d_air))

        if groove_depth < d_guide:
            etched = Slab(substrate(d_sub) + guiding(d_guide - groove_depth) \
			  + air(groove_depth + d_air))
        else:
            etched = Slab(substrate(d_sub - (groove_depth-d_guide)) \
			  + air(groove_depth + d_air))

        # Define 2D stack.

        stack = Stack(                                                       \
         waveguide(0) +                                                      \
         no_of_grooves * (      etched((1-filling_factor)* period_length)    \
                           + waveguide(   filling_factor * period_length)) + \
         waveguide(0))

        # Calculate.

        stack.calc()
    
        R = stack.R12(0,0)
        R_OK = 0.415900175747+0.0326513148086j
        print R, "expected", R_OK
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()
           
        self.failUnless(R_pass)

suite = unittest.makeSuite(grating2, 'test')

if __name__ == "__main__":
    unittest.main()
              
