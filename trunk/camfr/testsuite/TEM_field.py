#! /usr/bin/env python

###################################################################
#
# Test for field of TEM mode.
#
###################################################################

from camfr import *

import unittest, eps

class TEM_field(unittest.TestCase):
    def testTEM_field(self):
        
        """TEM field"""

        print
        print "Running TEM field..."

        set_N(1)
        set_polarisation(TM)

        set_lambda(1.55)

        GaAs = Material(3.5)
        air  = Material(1.0)

        s1 = Slab(air(0) + GaAs(4) + air(0))
        s2 = Slab(GaAs(2) + GaAs(2))
        s3 = Slab(GaAs(4))

        s1.calc()
        s2.calc()
        s3.calc()

        f1 = abs(s1.mode(0).field(Coord(3,0,0)).H2())
        f2 = abs(s2.mode(0).field(Coord(3,0,0)).H2())
        f3 = abs(s3.mode(0).field(Coord(3,0,0)).H2())
                
        f_OK = 0.0481935271885

        print f1, "expected", f_OK
        print f2, "expected", f_OK
        print f3, "expected", f_OK
                
        f1_pass = abs((f1 - f_OK) / f_OK) < eps.testing_eps
        f2_pass = abs((f2 - f_OK) / f_OK) < eps.testing_eps
        f3_pass = abs((f3 - f_OK) / f_OK) < eps.testing_eps
                
        free_tmps()

        set_polarisation(TE)
       
        self.failUnless(f1_pass and f2_pass and f3_pass)

suite = unittest.makeSuite(TEM_field, 'test')        

if __name__ == "__main__":
    unittest.main()



