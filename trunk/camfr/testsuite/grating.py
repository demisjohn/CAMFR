#! /usr/bin/env python

###################################################################
#
# Test for stable field calculation in a large grating.
#
###################################################################

from camfr import *
from Numeric import *
import unittest, eps

class grating(unittest.TestCase):
    def testgrating(self):

        """Grating"""

        print
        print "Running grating..."

        set_N(10)
        set_lambda(1.5)

        GaAs = Material(3.5)
        air  = Material(1.0)

        set_left_wall (slab_H_wall)
        set_right_wall(slab_H_wall)

        gp = 0.5*get_lambda()
        s1 = Slab(5*(air(gp/4.) + GaAs(gp/2.) + air(gp/4.))+air(0))
        s1.calc()
        E1 = s1.mode(0).field(Coord(gp/2.,0,0)).E2()

        E1_OK = -9.37774762123-2.02472587307e-14j
        print E1, "expected", E1_OK
        E1_pass = abs((E1 - E1_OK) / E1_OK) < eps.testing_eps


        gp = 1.0*get_lambda()
        s2 = Slab(5*(air(gp/4.) + GaAs(gp/2.) + air(gp/4.))+air(0))
        s2.calc()
        E2 = s2.mode(0).field(Coord(gp/2.,0,0)).E2()

        E2_OK = -7.03121528777-4.85905818077e-11j        
        print E2, "expected", E2_OK
        E2_pass = abs((E2 - E2_OK) / E2_OK) < eps.testing_eps

        free_tmps()
           
        self.failUnless(E1_pass and E2_pass)

suite = unittest.makeSuite(grating, 'test')

if __name__ == "__main__":
    unittest.main()
          
