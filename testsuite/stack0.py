#! /usr/bin/env python

####################################################################
#
# stack with zero length first slab
#
####################################################################

from camfr import *

import unittest, eps

class stack0(unittest.TestCase):
    def teststack0(self):
        
        """Stack0"""

        print('')
        print("Running stack0...")

        set_lambda(1.55)
        set_N(5)
        set_polarisation(TE)
        
        AlGaAs_m = Material(3.34)
        air_m    = Material(1)

        AlGaAs = Slab(AlGaAs_m(1))
        air    = Slab(air_m(1))
        
        s = Stack(AlGaAs(0) + air(1) + AlGaAs(0))

        inc = zeros(N())
        inc[0] = 1
        s.set_inc_field(inc)

        s.calc()
        
        E = s.field(Coord(0.5,0,0.5)).E2()
        E_OK = 37.3792383401-24.3729597789j
        
        print(E, "expected", E_OK)
        
        E_pass = abs((E - E_OK) / E_OK) < eps.testing_eps

        free_tmps()
        
        self.assertTrue(E_pass)

suite = unittest.makeSuite(stack0, 'test')        

if __name__ == "__main__":
    unittest.main()
    


