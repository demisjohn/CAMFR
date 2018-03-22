#! /usr/bin/env python

###################################################################
#
# Nearly degenerate rods.
#
###################################################################

from camfr import *

import unittest, eps

class rods(unittest.TestCase):
    def testrods(self):
        
        """Rods"""

        print('')
        print("Running rods...")
        
        set_lambda(1.5)
        set_N(100)
        set_mode_surplus(3)

        GaAs = Material(sqrt(11.4))
        air = Material(1.0)

        a = .600 # period
        r = .100 # rod radius

        set_upper_PML(-0.1)
        set_lower_PML(-0.1)
        
        d_clad = a

        periods = 5

        rods = Slab(  air(d_clad) + periods*(air(a-2*r) + GaAs(2*r)) \
                 + air(a) \
                 + periods*(air(a-2*r) + GaAs(2*r)) \
                 + air(a-2*r+d_clad) )

        rods.calc()

        n_eff_98    = rods.mode(98).n_eff()
        n_eff_98_OK = 0.520574469226-8.88879675497j

        print(n_eff_98, "expected", n_eff_98_OK)
        mode98_pass = abs((n_eff_98-n_eff_98_OK)/n_eff_98_OK) < eps.testing_eps

        n_eff_99    = rods.mode(99).n_eff()
        n_eff_99_OK = 0.520574469226-8.88879675497j

        print(n_eff_99, "expected", n_eff_99_OK)
        mode99_pass = abs((n_eff_99-n_eff_99_OK)/n_eff_99_OK) < eps.testing_eps
        
        free_tmps()

        set_upper_PML(0)
        set_lower_PML(0)

        self.failUnless(mode98_pass and mode99_pass)

suite = unittest.makeSuite(rods, 'test')        

if __name__ == "__main__":
    unittest.main()
    



