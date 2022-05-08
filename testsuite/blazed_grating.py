#!/usr/bin/env python

##############################################################################
#
# Blazed grating in slab waveguide. 
#
##############################################################################

from camfr import *
from cmath import *

import unittest, eps

class blazed_grating(unittest.TestCase):
    def testRT(self):
        
        """Blazed grating"""

        print('')
        print("Running blazed grating...")

        set_lambda(1.55)
        set_N(30)

        # Waveguide parameters.

        guide_thickness = .200
        groove_depth = 0.060
        d = 0.5
        
        set_lower_PML(-0.05)
        set_upper_PML(-0.05)
        
        # Grating parameters.

        no_of_grooves = 20
        period = 0.590
        filling_factor = 0.5

        # Blazing parameters

        theta = 0.5
        a = groove_depth * tan(theta)
        b = period/2.0 - a
        steps = 2

        # Materials.

        substrate = Material(1.5)
        guiding   = Material(3.5)
        air       = Material(1.0)

        # 1-D slab structures.

        waveguide = Slab(substrate(d) + guiding(guide_thickness) \
                         + air(d))

        etched = Slab(substrate(d) + \
                      guiding(guide_thickness - groove_depth) \
                      + air(groove_depth + d))

        # Left side of paralellogram.

        slices = []

        l_expr = Expression()

        for n in range(1,steps):
            delta_x = groove_depth*(steps-n)/steps
            slice = Slab(substrate(d)                           \
                 + guiding(guide_thickness - delta_x)           \
                 + air(delta_x + d) )
            slices.append(slice)
            l_expr.add(slice(a/steps))

        left = Stack(l_expr)

        # Right side of parallellogram.

        r_expr = Expression()

        for n in range(1,steps):
            slice = Slab(  substrate(d)                         \
                 + guiding(guide_thickness - groove_depth)      \
                 + air(groove_depth*n/steps)                    \
                 + guiding(groove_depth*(steps - n)/steps)      \
                 + air(d))      
            slices.append(slice)
            r_expr.add(slice(a/steps))

        right = Stack(r_expr)

        # Complete blazed grating.
        
        blazed = Stack(waveguide(period/2.0)                         \
             + no_of_grooves * ( Term(right) +    etched(b+a/steps)  \
                                +      left  + waveguide(b+a/steps)) \
             + waveguide(0))

        # Calc
        
        blazed.calc()

        R = pow(abs(blazed.R12(0,0)), 2)
        T = pow(abs(blazed.T12(0,0)), 2)

        R_OK = 0.400333359695
        T_OK = 0.204851452350

        print(R, "expected", R_OK)
        print(T, "expected", T_OK)

        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps
        T_pass = abs((T - T_OK) / T_OK) < eps.testing_eps

        free_tmps()

        set_lower_PML(0)
        set_upper_PML(0)

        self.assertTrue(R_pass and T_pass)

suite = unittest.makeSuite(blazed_grating, 'test')        

if __name__ == "__main__":
    unittest.main()
    
