#! /usr/bin/env python

###################################################################
#
# Calculates semi-infinite stack.
#
###################################################################

from camfr import *

from math import *

import unittest, eps

class blochstack(unittest.TestCase):
    def testblochstack(self):
        
        """Bloch stack"""

        print('')
        print("Running bloch stack...")

        set_N(40)
        set_polarisation(TE)

        GaAs = Material(3.37)
        air  = Material(1.00)

        a = 1.0
        r = 0.3 * a * sqrt(pi)/2.
        r_wg = 0.25 * a * sqrt(pi)/2.

        set_lower_wall(slab_H_wall)

        x_periods = 5

        set_upper_PML(-0.1)

        inc_wg_width = 0.5

        set_lambda(a/0.2)

        # Define Slabs.

        crystal = Slab(GaAs(r) + air(a-2*r) \
               + x_periods*(GaAs(2*r) + air(a-2*r)) + air(0))

        wg_1 = Slab(air(a-r) \
               + x_periods*(GaAs(2*r) + air(a-2*r)) + air(0))

        wg_2 = Slab(GaAs(r_wg) + air(a-r_wg-r) \
               + x_periods*(GaAs(2*r) + air(a-2*r)) + air(0))

        no_crystal = Slab(air(crystal.width()))

        inc_wg = Slab(GaAs(inc_wg_width/2.) \
               + air(crystal.width() - inc_wg_width/2.))

        # Calculate semi-infinite structure.

        period_0 = wg_1(r-r_wg) + wg_2(2*r_wg) + wg_1(r-r_wg) \
                   + no_crystal(a-2*r)

        s_inf = BlochStack(period_0)
        s = Stack(s_inf(0) + inc_wg(1) + s_inf(0))
        s.calc()

        R = s.R12(0,0)
        R_OK = -0.551728847141+0.66924073272j

        T = s.T12(0,1)
        T_OK = 0.0340145359674+0.154470138754j

        print(R, "expected", R_OK)
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        print(T, "expected", T_OK)
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        set_lower_wall(slab_E_wall)

        set_upper_PML(0)
       
        self.failUnless(R_pass and T_pass)

        free_tmps()

suite = unittest.makeSuite(blochstack, 'test')        

if __name__ == "__main__":
    unittest.main()
