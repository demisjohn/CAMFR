#! /usr/bin/env python

###################################################################
#
# Test for case where two waveguides are each other's mirror image..
#
###################################################################

from camfr_work import *

import unittest, eps

class shift(unittest.TestCase):
    def testsift(self):

        """Shift"""

        print
        print "Running shift..."

        set_N(2)
        set_polarisation(TE)
        set_lambda(1/0.7)

        rod = Material(sqrt(1))
        a = Material(sqrt(13.00))

        set_lower_wall(slab_H_wall)
        set_upper_wall(slab_H_wall)

        s1 = Slab(rod(0.43243175441)+a(0.0523712004609)+rod(0.381222448913))
        s2 = Slab(rod(0.381222448913)+a(0.0523712004609)+rod(0.43243175441))

        s = Stack(s1(1)+s2(1))
        s.calc()

        T = s.T12(0,0)
        T_OK = -0.141332685373-0.986658604736j
        print T, "expected", T_OK
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps
        
        free_tmps()
        
        set_upper_wall(slab_E_wall)
        set_lower_wall(slab_E_wall)
           
        self.failUnless(T_pass)

suite = unittest.makeSuite(shift, 'test')

if __name__ == "__main__":
    unittest.main()
          
