#! /usr/bin/env python

###################################################################
#
# Test for COST P11 photonic wire section..
#
###################################################################

from camfr_work import *

import unittest, eps

class section2(unittest.TestCase):
    def testsection2(self):

        """Section 2"""

        print
        print "Running section 2..."
        
        set_lambda(1.55)
        set_N(1)
        
        air_m  = Material(1.00)
        SiO2_m = Material(1.45)
        Si_m   = Material(3.50)

        set_left_wall(no_wall)
        set_right_wall(no_wall)

        W = 0.5 # width of wire
        D = 1.0 # height of ox buffer
        C = 1.5 # width of computational domain

        Si   = Slab(  Si_m(C))
        SiO2 = Slab(SiO2_m(C))
        air  = Slab( air_m(C))
        core = Slab(Si_m(W/2.)+air_m(C-W/2.))

        s = Section(core(.110)+SiO2(D)+Si(.1),core(.110)+air(.1), 50, 50)

        s.set_estimate(2.41237-2.9e-8j)
        s.calc()

        n_eff_0 = s.mode(0).n_eff()
        n_eff_0_OK = 2.41228943244-2.91399206483e-08j
        print n_eff_0, "expected", n_eff_0_OK
        n_eff_0_pass = abs((n_eff_0 - n_eff_0_OK)/n_eff_0_OK) < eps.testing_eps

        free_tmps()
        
        set_left_wall (E_wall)
        set_right_wall(E_wall)
        
        self.failUnless(n_eff_0_pass)

suite = unittest.makeSuite(section2, 'test')

if __name__ == "__main__":
    unittest.main()
