#! /usr/bin/env python

###################################################################
#
# Test various ways of creating expressions.
#
###################################################################

from camfr import *

import unittest

class expressions(unittest.TestCase):
    def test_expressions(self):
        
        """Expressions"""

        a = Material(1)
        waveguide = Slab(a(1))

        left = Stack(waveguide(1) + waveguide(2))

        wall = E_Wall(waveguide)
        e = Expression(waveguide(0) + left)

        e_tt = waveguide(0) + waveguide(1)
        e_te = waveguide(0) + e
        e_tw = waveguide(0) + wall
        e_ts = waveguide(0) + left

        e_et = e + waveguide(1)
        e_ee = e + e
        e_ew = e + wall
        e_es = e + left

        e_wt = wall + waveguide(1)
        e_we = wall + e
        e_ww = wall + wall
        e_ws = wall + left

        e_st = left + waveguide(1)
        e_se = left + e
        e_sw = left + wall
        e_ss = left + left

        self.failUnless(1)

suite = unittest.makeSuite(expressions, 'test')        

if __name__ == "__main__":
    unittest.main()
        
