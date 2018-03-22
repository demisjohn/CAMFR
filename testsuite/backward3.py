#! /usr/bin/env python

############################################################################
#
# Backward modes 3.
#
#  Note: this one has backward modes, but switching them to forward modes
#  introduces gain. We pick the sign for lossy propagation
#
############################################################################

from camfr import *

import unittest, eps

class backward3(unittest.TestCase):
    def testbackward3(self):

        """backward 3"""

        print('')
        print("Running backward 3...")

        set_lambda(1.64367524538)
        set_N(20)
        set_solver(stretched_ASR)
        set_orthogonal(False)
        set_mode_surplus(10)
        set_polarisation(TM)
        set_upper_wall(slab_H_wall)
        set_lower_wall(slab_H_wall)

        mat_hi = Material(2.0)
        mat_lo = Material(1.0)
        metal  = Material(-1.43028859905j)

        wavelen_plasma = 0.941825915601
        
        d       = 0.02*wavelen_plasma 
        d_clad  = 1.70*wavelen_plasma
        d_metal = 1.20*wavelen_plasma
    
        slab1 = Slab(metal(d_metal)+mat_hi(d)+mat_lo(d_clad-d))
        slab2 = Slab(metal(d_metal)+mat_hi(d_clad))
        stack = Stack(slab1(0)+slab2(0))

        stack.calc()
    
        R = stack.R12(0,0)
        R_OK = -0.246126984037-0.00532954414873j
        print(R, "expected", R_OK)
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps
        
        set_solver(track)
        set_mode_surplus(1.5)
    
        set_polarisation(TE)
        set_upper_wall(slab_E_wall)
        set_lower_wall(slab_E_wall)
        set_orthogonal(True)
        
        self.failUnless(R_pass)

suite = unittest.makeSuite(backward3, 'test')

if __name__ == "__main__":
    unittest.main()
