#! /usr/bin/env python

###################################################################
#
# Test for coupled waveguides.
#
###################################################################

from camfr import *

import unittest, eps

class coupled(unittest.TestCase):
    def testcoupled(self):

        """Coupled"""

        print('')
        print("Running coupled...")

        set_N(4)
        set_mode_surplus(10)
        
        set_lambda(1.55)
        set_polarisation(TE)

        air     = Material(1.0)
        BCB     = Material(1.53)
        SiO2    = Material(1.45)
        InP_eff = Material(3.1)
        Si_eff  = Material(2.728)

        s = Slab(SiO2(2.0)+Si_eff(0.22)+BCB(.5)+InP_eff(0.15)+air(2.0))
        s.calc()

        n_eff_0 = s.mode(0).n_eff()
        n_eff_0_OK = 2.1914030141
        print(n_eff_0, "expected", n_eff_0_OK)
        n_eff_0_pass = abs((n_eff_0 - n_eff_0_OK)/n_eff_0_OK) < eps.testing_eps

        n_eff_1 = s.mode(1).n_eff()
        n_eff_1_OK = 2.14181217076
        print(n_eff_1, "expected", n_eff_1_OK)
        n_eff_1_pass = abs((n_eff_1 - n_eff_1_OK)/n_eff_1_OK) < eps.testing_eps

        free_tmps()
           
        self.assertTrue(n_eff_0_pass and n_eff_1_pass)

suite = unittest.makeSuite(coupled, 'test')

if __name__ == "__main__":
    unittest.main()
