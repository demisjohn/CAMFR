#! /usr/bin/env python

###################################################################
#
# Test for COST P11 photonic wire section, TM mode.
#
###################################################################

from camfr_work import *

import unittest, eps

class section3(unittest.TestCase):
    def testsection3(self):

        """Section 3"""

        print
        print "Running section 3..."
        
        set_lambda(1.55)
        set_N(1)
        
        air_m  = Material(1.00)
        SiO2_m = Material(1.45)
        Si_m   = Material(3.50)

        set_left_wall(no_wall)
        set_right_wall(no_wall)
        
        set_upper_wall(slab_H_wall)
        set_lower_wall(slab_H_wall)
        
        W = 0.5 # width of wire
        D = 1.0 # height of ox buffer
        C = 1.5 # width of computational domain

        Si   = Slab(  Si_m(C))
        SiO2 = Slab(SiO2_m(C))
        air  = Slab( air_m(C))
        core = Slab(Si_m(W/2.)+air_m(C-W/2.))

        s = Section(core(.110)+SiO2(D)+Si(.1),core(.110)+air(.1), 50, 50)

        s.set_estimate(1.6-0.00066j)
        s.calc()

        n_eff_0 = s.mode(0).n_eff()
        n_eff_0_OK = 1.5999984535-0.000665487870832j
        print n_eff_0, "expected", n_eff_0_OK
        n_eff_0_pass = abs((n_eff_0 - n_eff_0_OK)/n_eff_0_OK) < eps.testing_eps

        free_tmps()
        
        set_left_wall (E_wall)
        set_right_wall(E_wall)
        
        set_upper_wall(slab_E_wall)
        set_lower_wall(slab_E_wall)
    
        self.failUnless(n_eff_0_pass)

suite = unittest.makeSuite(section3, 'test')

if __name__ == "__main__":
    unittest.main()
