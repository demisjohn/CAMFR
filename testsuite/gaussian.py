#! /usr/bin/env python

###################################################################
#
# Gaussian incident field.
#
###################################################################

from camfr import *

import unittest, eps

class gaussian(unittest.TestCase):
    def testgaussian(self):

        """Gaussian"""

        print('')
        print("Running Gaussian...")

        set_N(40)
        set_lambda(1.55)
        set_mode_surplus(2)
        
        GaAs = Material(3.5)
        air  = Material(1.0)

        set_lower_PML(-0.1)
        set_upper_PML(-0.1)
        
        slab1 = Slab(air(2)+ GaAs(1)+ air(2))
        slab2 = Slab(air(slab1.width()))
    
        s = Stack(slab1(0)+slab2(0))
        s.set_inc_field_gaussian(4,0.5,2.5, 0.001)

        s.calc()

        T = abs(s.trans_field()[0])
        T_OK = 0.217622614347

        print(T, "expected", T_OK)
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        free_tmps()

        set_lower_PML(0)
        set_upper_PML(0)
           
        self.failUnless(T_pass)

suite = unittest.makeSuite(gaussian, 'test')

if __name__ == "__main__":
    unittest.main()
