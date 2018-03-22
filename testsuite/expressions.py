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

        print('')
        print("Running expressions...")

        a = Material(1)
        
        s  = Slab(a(1))
        s2 = Slab(a(1))
        
        st = Stack(s(1) + s(2))

        wall = E_Wall(s)
        e = Expression(s(0) + st)

        e_tt = s(0) + s(1)
        e_te = s(0) + e
        e_tw = s(0) + wall
        e_ts = s(0) + st

        e_et = e + s(1)
        e_ee = e + e
        e_ew = e + wall
        e_es = e + st

        e_wt = wall + s(1)
        e_we = wall + e
        e_ww = wall + wall
        e_ws = wall + st

        e_st = st + s(1)
        e_se = st + e
        e_sw = st + wall
        e_ss = st + st

        stack = Stack(s(0))
        stack = Stack(s(1) + s(1))
        stack = Stack(s(0) + s(0))
        stack = Stack(s(0) + s(1))
        stack = Stack(s(0) + s(1) + s(0))
        stack = Stack(s(0) + s(0) + s(0))

        stack = Stack(s2(0) + s(0))
        stack = Stack(s2(0) + s(1) + s(1))
        stack = Stack(s2(0) + s(0) + s(0))
        stack = Stack(s2(0) + s(0) + s(1))
        stack = Stack(s2(0) + s(0) + s(1) + s(0))
        stack = Stack(s2(0) + s(0) + s(0) + s(0))

        stack = Stack(s2(1) + s(0))
        stack = Stack(s2(1) + s(1) + s(1))
        stack = Stack(s2(1) + s(0) + s(0))
        stack = Stack(s2(1) + s(0) + s(1))
        stack = Stack(s2(1) + s(0) + s(1) + s(0))
        stack = Stack(s2(1) + s(0) + s(0) + s(0))

        free_tmps()

        self.failUnless(1)

suite = unittest.makeSuite(expressions, 'test')        

if __name__ == "__main__":
    unittest.main()
        
