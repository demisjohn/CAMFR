#! /usr/bin/env python

##############################################################################
#
# Test various complicated ways of forming expressions.
#
##############################################################################

from camfr_work import *

import unittest, eps

class substacks(unittest.TestCase):
    def testRT(self):
        
        """Substracks"""

        print
        print "Running substacks..."

        set_lambda(1.1)
        set_N(50)

        a = Material(1)
        b = Material(2)

        e1 = Expression(a(1.0) + b(0.5))
        e2 = Expression(a(0.5) + b(1.0))

        s0 = Slab(a(1) + e1 + e2 + a(1))

        s1 = Slab(e1 + 2*e2 + a(0.5))

        e3 = Expression(s0(1) + s1(2))

        st = Stack(e3 + s0(1))

        st.calc()

        R = abs(st.R12(0,0))

        R_OK = 0.155960743972

        print R, "expected", R_OK

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(R_pass)


suite = unittest.makeSuite(substacks, 'test')        


if __name__ == "__main__":
    unittest.main()
