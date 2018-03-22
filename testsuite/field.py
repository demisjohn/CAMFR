#! /usr/bin/env python

####################################################################
#
# field profiles
#
####################################################################

from camfr import *

import unittest, eps

class field(unittest.TestCase):
    def testfield(self):
        
        """Field"""

        print('')
        print("Running field...")

        set_lambda(1.0)
        set_N(40)
        PML = -0.1
        set_lower_PML(PML)
        set_upper_PML(PML)
        set_polarisation(TE)
 
        wgcore = Material(3.5)
        wgclad = Material(2.5)
        wg = Slab(wgclad(2.25) + wgcore(0.5) + wgclad(2.25))
 
        stack=Stack(wg(2.0))
        stack.calc()

        f1 = wg.mode(1).field(Coord(0,0,0)).E2()
        f_OK = 5.30053306821e-11-2.14419512113e-11j

        print(f1, "expected", f_OK)
        f1_pass = abs(f1 - f_OK) < eps.testing_eps
        
        inc = zeros(N())
        inc[1] = 1
        stack.set_inc_field(inc)

        f2 = stack.field(Coord(0,0,0)).E2()

        print(f2, "expected", f_OK)
        f2_pass = abs(f2 - f_OK) < eps.testing_eps       

        free_tmps()

        set_upper_PML(0)
        set_lower_PML(0)
        
        self.failUnless(f1_pass and f2_pass)

suite = unittest.makeSuite(field, 'test')        

if __name__ == "__main__":
    unittest.main()
