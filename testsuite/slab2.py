#! /usr/bin/env python

###################################################################
#
# Slab.
#
###################################################################

from camfr import *

import unittest, eps

class slab2(unittest.TestCase):
    def testslab2(self):

        """Slab 2"""

        print('')
        print("Running slab 2...")
        
        set_N(1)

        set_upper_PML(-0.05)
        set_lower_wall(slab_H_wall)

        set_lambda(1./.2)

        GaAs = Material(3.37)
        air  = Material(1.0)

        s = Slab(GaAs(0.159358)+air(5.84064))

        s.calc()

        f = s.mode(0).field(Coord(0,0,0)).E2()
        f_OK = 24.6737+1.24772e-09j

        print(f, "expected", f_OK)
        f_pass = abs((f - f_OK) / f_OK) < eps.testing_eps

        free_tmps()

        set_upper_PML(0)

        set_lower_wall(slab_E_wall)
           
        self.failUnless(f_pass)

suite = unittest.makeSuite(slab2, 'test')

if __name__ == "__main__":
    unittest.main()
