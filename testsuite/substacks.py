#! /usr/bin/env python

##############################################################################
#
# Test various complicated ways of forming expressions.
#
##############################################################################

from camfr import *

import unittest, eps

class substacks(unittest.TestCase):
    def testRT(self):
        
        """Substacks"""

        print('')
        print("Running substacks...")

        set_lambda(1.1)
        set_N(50)

        set_lower_PML(-0.001)
        set_upper_PML(-0.001)

        a = Material(1)
        b = Material(2)

        e1 = Expression(a(0.6) + b(0.5))
        e2 = Expression(a(0.5) + b(0.6))

        s0 = Slab(a(0.2) + e1 + e2 + a(1.2))

        s1 = Slab(a(0.1) + e1 + 2*e2 + a(0.2))

        e3 = Expression(s0(1) + s1(2))

        st = Stack(e3 + s0(1))
        
        st.calc()

        R = abs(st.R12(0,0))

        R_OK = 0.0362401973292

        print(R, "expected", R_OK)

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(R_pass)


suite = unittest.makeSuite(substacks, 'test')        


if __name__ == "__main__":
    unittest.main()
