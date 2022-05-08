#! /usr/bin/env python

#############################################################################
#
# Polaritonic material.
#
#############################################################################

from camfr import *
from cmath import *

import unittest, eps

class polariton2(unittest.TestCase):
    def testpolariton2(self):

        """Polariton 2"""

        print('')
        print("Running polariton 2...")

        set_N(50)
        set_polarisation(TM)
        set_solver(series)
        set_mode_surplus(4)
        set_keep_all_1D_estimates(0)

        # Define structure

        z_periods = .5

        a = 1.0
        r = 0.25 * a

        set_lambda(a/.44)

        air = Material(1.0)
        polariton = Material(1.0999937705-10.8034009584j)
        
        g = Geometry(air)
        g += Square(Point(0,0.5*a), r, polariton)

        # Create expression.

        dx = a/40.
        dy = a/40.

        e = g.to_expression(-a,  a, dx, 0,  a, dy)

        s = Stack(e)
        s.calc()

        R = s.R12(0,0)
        R_OK = 0.129865230383-0.18515631071j
        
        print(R, "expected", R_OK)
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()
        
        set_polarisation(TE)
        set_solver(track)
        set_mode_surplus(1.5)
        set_keep_all_1D_estimates(0)
        
        self.assertTrue(R_pass)

suite = unittest.makeSuite(polariton2, 'test')

if __name__ == "__main__":
    unittest.main()
