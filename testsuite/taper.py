#! /usr/bin/env python

###################################################################
#
# Calculates a taper taper created from a geometry.
#
###################################################################

from camfr import *
from math import *

import unittest, eps

class taper(unittest.TestCase):
    def testtaper(self):
        
        """Taper"""

        print
        print "Running taper..."
        
        set_N(40)

        GaAs = Material(3.37)
        air  = Material(1.0)

        # Define structure

        a = 1.0
        r = 0.3 * a * sqrt(pi)/2.
        r_wg = 0.25 * a * sqrt(pi)/2.

        set_lower_wall(slab_H_wall)

        set_upper_PML(-0.1)

        wg_width = 1 # in units of a
        taper_length = 4 # in units of a

        x_periods = 4
        z_periods = taper_length + 2

        # Create geometry.

        g = Geometry(air)

        for ix in range(0,z_periods+1):
          g += Square(Point(ix*a,0), 2*r_wg, GaAs)
          for iy in range(1,x_periods+1):
            g += Square(Point(ix*a,iy*a), 2*r, GaAs)
    
        g += Rectangle(Point(-2*r-2*a,0), Point(-r, wg_width*a/2.), GaAs)
        g += Triangle(Point(-r,0), Point(-r, wg_width*a/2.),
                      Point(taper_length*a-r_wg,r_wg),GaAs)
        g += Triangle(Point(-r,0), Point(taper_length*a-r_wg,0),
                      Point(taper_length*a-r_wg,r_wg),GaAs) 

        dx = a/20.
        dy = a/20.

        e = Expression(g.to_expression(-2*r-2*a,     z_periods*a, dx,
                                           0,    (x_periods+1)*a, dy))

        s = Stack(e)

        e_inf = Expression(g.to_expression(z_periods*a - a, z_periods*a, dx,\
                                                    0,  (x_periods+1)*a, dy))
            
        s_inf = InfStack(e_inf)

        rx = arange(0.0, (x_periods+1)*a, a/10.)
        rz = arange(0.0, s.length(),      a/10.)

        taper = Stack(e + s_inf)

        # Calculate reflection.
            
        set_lambda(a/.2335)

        taper.calc()

        R = taper.R12(0,0)
        R_OK = 0.0273696979612+0.337022950102j
    
        print R, "expected", R_OK
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        set_lower_wall(slab_E_wall)

        set_upper_PML(0)
       
        self.failUnless(R_pass)

suite = unittest.makeSuite(taper, 'test')        

if __name__ == "__main__":
    unittest.main()

