#! /usr/bin/env python

###################################################################
#
# Testcase for ADR complex root finder.
#
###################################################################

from camfr_work import *

import unittest, eps

class ADR_solver(unittest.TestCase):
    def testADRsolver(self):
        
        """ADR solver"""

        print
        print "Running ADR solver..."
        
        set_solver(ADR)
        set_N(10)
        set_polarisation(TM)
        set_lambda(1.55)

        FeCo         = Material(3.565458-7.13499488050650j)
        InP          = Material(3.18)
        InGaAs       = Material(3.6-0.0806676j)
        InGaAsP_1_25 = Material(3.36)
        InGaAsP_1_55 = Material(3.61)

        s_ADR = Slab(InP(1.5) + InGaAsP_1_25(0.1) + InGaAsP_1_55(.15)    \
		     + InGaAsP_1_25(0.1)                                 \
		     + InP(0.5) + InGaAs(0.05) + FeCo(0.05))

        s_ADR.calc()

        n_eff_0    = s_ADR.mode(0).n_eff()
        n_eff_0_OK = 3.51054-0.349863j

        print n_eff_0, "expected", n_eff_0_OK
        mode0_pass = abs((n_eff_0 - n_eff_0_OK)/n_eff_0_OK) < eps.testing_eps
        
        free_tmps()
        set_solver(track)

        self.failUnless(mode0_pass)

suite = unittest.makeSuite(ADR_solver, 'test')        

if __name__ == "__main__":
    unittest.main()
    

