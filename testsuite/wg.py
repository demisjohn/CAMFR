#!/usr/bin/env python

####################################################################
#
# Sign of kz.
#
####################################################################

from camfr import *

import unittest, eps

class wg(unittest.TestCase):
    def testwg(self):
        
        """Waveguide"""

        print
        print "Running waveguide..."

        set_lambda(1)
        set_N(20)

        air = Material(1)
        mat = Material(3)

        wg = Slab(air(2-.1j)+mat(1)+air(2-.1j))
        wg.calc()

        n =  wg.mode(0).n_eff()
        n_OK = 2.96617486466+9.44438326765e-17j

        print n, "expected", n_OK
        n_pass = abs((n - n_OK) / n_OK) < eps.testing_eps

        free_tmps()
       
        self.failUnless(n_pass)

suite = unittest.makeSuite(wg, 'test')        

if __name__ == "__main__":
    unittest.main()
    
     
