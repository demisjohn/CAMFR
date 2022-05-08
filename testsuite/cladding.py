#! /usr/bin/env python

###################################################################
#
# Test for stable field calculation in cladding.
#
###################################################################

from camfr import *

import unittest, eps

class cladding(unittest.TestCase):
    def testcladding(self):

        """cladding"""

        print('')
        print("Running cladding...")
        
        set_N(4)
        set_polarisation(TE)
        set_mode_surplus(40)

        set_lambda(1.5)

        GaAs = Material(3.5)
        air  = Material(1.0)
      
        core = Slab(air(10.1)+ GaAs(1)+ air(10.1))
        core.calc()
      
        E_field = core.mode(0).field(Coord(core.width()-0.01,0,0)).E2()
        E_field_OK = 0.0
      
        print(E_field, "expected", E_field_OK)

        E_field_pass = abs(E_field) < eps.testing_eps

        free_tmps()
      
        self.assertTrue( E_field_pass )

suite = unittest.makeSuite(cladding, 'test')        

if __name__ == "__main__":
    unittest.main()
          
