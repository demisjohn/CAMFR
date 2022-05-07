#! /usr/bin/env python

###################################################################
#
# Backward modes 2.
#
###################################################################

from camfr import *

import unittest, eps

class backward2(unittest.TestCase):
    def testbackward2(self):

        """backward 2"""

        print('')
        print("Running backward 2...")

        a = 0.5
    
        set_lambda(a/0.183)  
        set_N(10) 
        set_backward_modes(1)
        set_circ_PML(-0.0001)
        
        air = Material(1.0)
        m1  = Material(3.5)
        m2  = Material(1.45)
    
        r      = 0.62*a
        d_clad = a-r
        d1     = 0.4*a
        d2     = 0.95*a
    
        c = Circ(m1(r)+air(d_clad)+m1(d1)+m2(d2)+m1(d1)+m2(d2)+m1(d1)+m2(d2))
        c.calc()
    
        n_eff = c.mode(7).n_eff()
        n_eff_OK = -1.45872127197e-05-0.123469672815j
        print(n_eff, "expected", n_eff_OK)
        n_eff_pass = abs((n_eff - n_eff_OK) / n_eff_OK) < eps.testing_eps
        
        set_circ_PML(0)   
        set_backward_modes(0)
           
        self.assertTrue(n_eff_pass)

suite = unittest.makeSuite(backward2, 'test')

if __name__ == "__main__":
    unittest.main()    
