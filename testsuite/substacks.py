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
        
        """Substracks"""

        set_lambda(1.1)
        set_N(50)

        a = Material(1)
        b = Material(2)

        e1 = Expression(a(2)+b(1))
        e2 = Expression(a(1)+b(2))

        e1 + Term(e2)

        s0 = Slab(a(2) + e1 + Term(e2) + a(1))

        s1 = Slab(a(2) + Term(e1) + Term(e2))

        s2 = Slab(e1+2*Term(e2))

        e3 = Expression(s0(1) + s2(2))

        st = Stack(e3 + s0(1))

        st.calc()

        R = abs(st.R12(0,0))

        R_OK = 0.63344859662

        print R, "expected", R_OK

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(R_pass)


suite = unittest.makeSuite(substacks, 'test')        


if __name__ == "__main__":
    unittest.main()
