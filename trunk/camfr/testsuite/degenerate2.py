#! /usr/bin/env python

###################################################################
#
# Test for uncoupled waveguides.
#
###################################################################

from camfr import *

import unittest, eps

class degenerate2(unittest.TestCase):
    def testdegenerate2(self):

        """Degenerate 2"""

        print
        print "Running degenerate 2..."

        set_lambda(1.55)
        set_N(8)

        set_field_calc_heuristic(symmetric)
        set_unstable_exp_threshold(1e-6)

        a = Material(1.0)
        M = Material(3.4)
        m = Material(3.08221)

        s1 = Slab(a(0.20571)+M(0.15995)+a(0.41142)+M(0.15995)+a(0.98279)+\
           m(0.15995)+a(0.41142)+M(0.15995)+a(0.98279)+M(0.15995)+a(0.41142)+\
           m(0.15995)+a(0.98279)+M(0.15995)+a(0.41142)+M(0.15995)+a(0.20571))

        s2 = Slab(a(s1.width()))

        s = Stack(s2(0)+s1(0))
        s.calc()

        T = s.T12(0,0)
        T_OK = 0.275422970158+9.09323340151e-05j

        print T, "expected", T_OK
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        free_tmps()

        set_unstable_exp_threshold(1e-12)
        set_field_calc_heuristic(identical)
        
        self.failUnless(T_pass)

suite = unittest.makeSuite(degenerate2, 'test')

if __name__ == "__main__":
    unittest.main()
          
