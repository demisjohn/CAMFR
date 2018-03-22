#! /usr/bin/env python

###################################################################
#
# fw_bw.
#
###################################################################

from camfr import *

import unittest, eps

class fw_bw(unittest.TestCase):
    def testfw_bw(self):

        """fw_bw"""

        print('')
        print("Running fw_bw...")
	
        set_lambda(1)
        set_N(20)
        set_polarisation(TE)
        set_lower_PML(-0.1)
        set_upper_PML(-0.1)
        
        core = Material(3.0)
        clad = Material(1.0)

        s1 = Slab(clad(1) + core(.5) + clad(1))
        s2 = Slab(clad(s1.width()))

        s = Stack(s1(0)+s2(1)+s1(0))

        inc = zeros(N())
        inc[0] = 1
        s.set_inc_field(inc)
        
        s.calc()

        fw, bw = s.fw_bw(.5)

        f = abs(fw[0])
        f_OK = abs(0.556349919991+0.0468263510045j)
        print(f, "expected", f_OK)
        f_pass = abs((f - f_OK) / f_OK) < eps.testing_eps
        
        b = abs(bw[0])
        b_OK = abs(-0.139412386777-0.0500728037444j)
        print(b, "expected", b_OK)
        b_pass = abs((b - b_OK) / b_OK) < eps.testing_eps

        free_tmps()

        set_upper_PML(0)
        set_lower_PML(0)
           
        self.failUnless(f_pass and b_pass)

suite = unittest.makeSuite(fw_bw, 'test')

if __name__ == "__main__":
    unittest.main()
