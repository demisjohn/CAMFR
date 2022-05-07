#! /usr/bin/env python

###################################################################
#
# Test for coupled waveguides in a section geometry.
#
###################################################################

from camfr import *

import unittest, eps

class section1(unittest.TestCase):
    def testsection1(self):

        """Section 1"""

        print('')
        print("Running section 1...")
        
        set_lambda(1.55)
        set_N(2)

        core = Material(3.47572)
        clad = Material(1.44402)
        air = Material(1)

        set_lower_wall(slab_H_wall)
        set_upper_wall(slab_H_wall)

        d   = 0.22
        w   = 0.43
        gap = 0.18
        
        wg1 = Slab(clad(1) + core(d) + air(1))
        wg2 = Slab(clad(1) + air(d+1))

        s = Section(wg2(1)+wg1(w)+wg2(gap)+wg1(w)+wg2(1), 10, 40)
        
        s.set_estimate(2.20)
        s.set_estimate(1.59)
        s.calc()

        n_eff_0 = s.mode(0).n_eff()
        n_eff_0_OK = 2.20461590777
        print(n_eff_0, "expected", n_eff_0_OK)
        n_eff_0_pass = abs((n_eff_0 - n_eff_0_OK)/n_eff_0_OK) < eps.testing_eps

        n_eff_1 = s.mode(1).n_eff()
        n_eff_1_OK = 1.59844475105
        print(n_eff_1, "expected", n_eff_1_OK)
        n_eff_1_pass = abs((n_eff_1 - n_eff_1_OK)/n_eff_1_OK) < eps.testing_eps

        free_tmps()
        
        set_lower_wall(slab_E_wall)
        set_upper_wall(slab_E_wall)
        
        self.assertTrue(n_eff_0_pass and n_eff_1_pass)

suite = unittest.makeSuite(section1, 'test')

if __name__ == "__main__":
    unittest.main()
