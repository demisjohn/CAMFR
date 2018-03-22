#! /usr/bin/env python

###################################################################
#
# Backward modes.
#
###################################################################

from camfr import *

import unittest, eps

class backward(unittest.TestCase):
    def testbackward(self):

        """backward"""

        print('')
        print("Running backward...")


        mat = Material(3.5)
        air = Material(1.0)

        wvg = Circ(mat(0.6) + air(0.4))
        sub = Circ(mat(1))

        set_N(20)
        set_lambda(1/0.18)
        set_circ_order(1)
        set_backward_modes(1)
        
        s = Stack(wvg(1)+sub(1))
        s.calc()

        R = s.R12(0,0)
        R_OK = -0.0392923220796+0.0408718742985j
        print(R, "expected", R_OK)
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        T = s.T12(0,0)
	T_OK = 0.202336029811+0.776634435067j
        print(T, "expected", T_OK)
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        set_backward_modes(0)
           
        self.failUnless(R_pass and T_pass)

suite = unittest.makeSuite(backward, 'test')

if __name__ == "__main__":
    unittest.main()
