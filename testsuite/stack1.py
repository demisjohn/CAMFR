#! /usr/bin/env python

####################################################################
#
# stack with zero length first slab
#
####################################################################

from camfr import *

import unittest, eps

class stack1(unittest.TestCase):
    def teststack1(self):
        
        """Stack1"""

        print('')
        print("Running stack1...")

        set_N(40)
        set_lambda(1.55)

        # Define materials.

        GaAs = Material(3.5)
        air  = Material(1.0)

        # Define geometry.

        set_lower_PML(-0.1)
        set_upper_PML(-0.1)

        wg  = Slab(air(2.0) + GaAs(0.2) + air(2.0))
        gap = Slab(air(4.2))

        s = Stack(wg(0) + gap(1) + wg(0))

        inc = zeros(N())
        inc[0] = 1
        s.set_inc_field(inc)

        n = s.n(Coord(2.1,0,-1))
        n_OK = 3.5
        
        print(n, "expected", n_OK)

        n_pass = abs((n - n_OK) / n_OK) < eps.testing_eps

        E = abs(s.field(Coord(2.1,0,-1)).E2())
        E_OK = 11.0402395717
        
        print(E, "expected", E_OK)
        
        E_pass = abs((E - E_OK) / E_OK) < eps.testing_eps

        free_tmps()
        
        self.failUnless(n_pass and E_pass)

suite = unittest.makeSuite(stack1, 'test')        

if __name__ == "__main__":
    unittest.main()
