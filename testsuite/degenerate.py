#! /usr/bin/env python

###################################################################
#
# Test for uncoupled waveguides.
#
###################################################################

from camfr import *

import unittest, eps

class degenerate(unittest.TestCase):
    def testdegenerate(self):

        """Degenerate"""

        print('')
        print("Running degenerate...")

        set_N(20)
        set_lambda(1.00)
        set_polarisation(TE)
        set_mode_surplus(7)

        m = Material(3.5)
        a = Material(1.0)

        s = Slab(a(2) + m(0.5) + a(2.99) + m(0.5) + a(2))
        s.calc()
        
        E0_l = abs(s.mode(0).field(Coord(2.25,0,0)).E2())
        E0_r = abs(s.mode(0).field(Coord(s.width()-2.25,0,0)).E2())
        
        E1_l = abs(s.mode(1).field(Coord(2.25,0,0)).E2())
        E1_r = abs(s.mode(1).field(Coord(s.width()-2.25,0,0)).E2())
        
        E_OK = 19.2559587762

        OK = 0

        if abs(E0_l) > abs(E0_r): # Mode 0 in left waveguide
            OK =   (abs((E0_l - E_OK) / E_OK) < eps.testing_eps) \
               and (abs((E1_r - E_OK) / E_OK) < eps.testing_eps) \
               and (abs(E0_r) < eps.testing_eps) \
               and (abs(E1_l) < eps.testing_eps)
        else: # Mode 1 in left waveguide
            OK =   (abs((E1_l - E_OK) / E_OK) < eps.testing_eps) \
               and (abs((E0_r - E_OK) / E_OK) < eps.testing_eps) \
               and (abs(E1_r) < eps.testing_eps) \
               and (abs(E0_l) < eps.testing_eps)

        free_tmps()
           
        self.assertTrue(OK)

suite = unittest.makeSuite(degenerate, 'test')

if __name__ == "__main__":
    unittest.main()
          
