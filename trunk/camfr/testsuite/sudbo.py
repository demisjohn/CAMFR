#!/usr/bin/env python

####################################################################
#
# Reproduces 2D waveguide simulations from Sudbo's paper.
#
####################################################################

from camfr import *

import unittest, eps

class sudbo(unittest.TestCase):
    def testsudbo(self):
        
        """Sudbo"""

        print
        print "Running Sudbo..."

        set_lambda(1.0)
        set_N(1)

        core = Material(1.5)
        clad = Material(1.0)

        set_mode_correction(full)

        set_left_wall(H_wall)
        set_right_wall(H_wall)

        wg = Slab(core(0.25) + clad(0.76-.25))
        air = Slab(clad(wg.width()))
        
        all = wg(0.25) + air(3.01-.25)

        s = Section(all, 6, 40)
        s.calc()

        n_eff = s.mode(0).n_eff()
        n_eff_OK = 1.16243992419
        
        print n_eff, "expected", n_eff_OK
        n_eff_pass = abs((n_eff - n_eff_OK)/n_eff_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(n_eff_pass)

suite = unittest.makeSuite(sudbo, 'test')        

if __name__ == "__main__":
    unittest.main()
    
