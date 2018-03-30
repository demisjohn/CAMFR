#! /usr/bin/env python

###################################################################
#
# Slab.
#
###################################################################

from camfr import *

import unittest, eps

class slab(unittest.TestCase):
    def testslab(self):

        """Slab"""

        print('')
        print("Running slab...")
    
        set_lambda(1.56)
        set_N(75)
        set_mode_surplus(2)
        set_polarisation(TM)
        set_upper_PML(-0.25)   
        core = Material(2.8305)
        clad = Material(1.0)

        #set_mode_surplus(3)

        s = Slab(core(4.3) + clad(12.95))
        s.calc()

        n = s.mode(74).n_eff()
        n_OK = 0.0130942-2.94038j
    
        print(n, "expected", n_OK)
        n_pass = abs((n - n_OK) / n_OK) < eps.testing_eps

        free_tmps()

        set_upper_PML(0)
        set_lower_PML(0)
        set_polarisation(TE)
           
        self.failUnless(n_pass)

suite = unittest.makeSuite(slab, 'test')

if __name__ == "__main__":
    unittest.main()
