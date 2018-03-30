#! /usr/bin/env python

###################################################################
#
# Slab.
#
###################################################################

from camfr import *

import unittest, eps

class slab3(unittest.TestCase):
    def testslab3(self):

        """Slab 3"""

        print('')
        print("Running slab 3...")
        
        set_N(10)

        set_upper_PML(-0.08)
        set_lower_wall(slab_H_wall)

        set_lambda(1.545)

        m1 = Material(3.33211366381)
        m2 = Material(3.22559213858507)
        m3 = Material(1.49097599017)
        s = Slab(m1(0.143333333333)+m2(0.02866666666)+m3(2.838))

        s.calc()

        n0 = s.mode(0).n_eff()
        n0_OK = 2.96986800836

        print(n0, "expected", n0_OK)
        n0_pass = abs((n0 - n0_OK) / n0_OK) < eps.testing_eps

        n9 = s.mode(9).n_eff()        
        n9_OK = 0.0814435888385-1.84709216839j
 
        print(n9, "expected", n9_OK)
        n9_pass = abs((n9 - n9_OK) / n9_OK) < eps.testing_eps

        free_tmps()

        set_upper_PML(0)

        set_lower_wall(slab_E_wall)
           
        self.failUnless(n0_pass and n9_pass)

suite = unittest.makeSuite(slab3, 'test')

if __name__ == "__main__":
    unittest.main()
