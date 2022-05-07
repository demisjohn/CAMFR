#! /usr/bin/env python

###################################################################
#
# Calculates modification of spontaneous emission of a dipole
# between two metal plates.
#
###################################################################

from camfr import *
from math import *

import unittest, eps

class SpE(unittest.TestCase):
    def testSpE(self):
        
        """Spontaneous emission"""

        print('')
        print("Running spontaneous emission...")

        set_N(60)
        set_lambda(1)
        set_circ_order(1)
        set_circ_field_type(cos_type)
        set_circ_PML(-0.5)

        # Define waveguide and wall.

        air_m = Material(1.0)
        air = Circ(air_m(10))

        wall = E_Wall(air)

        # Define cavities.

        d = 2.75

        s      = Stack(air(d/2.) + wall)
        s_open = Stack(air(d))

        source_pos  = Coord(0,0,0)
        orientation = Coord(1,0,0)

        # Calculate modification of spontaneous emission.
    
        cav = Cavity(s, s)
        cav.set_source(source_pos, orientation)

        cav_open = Cavity(s_open, s_open)
        cav_open.set_source(source_pos, orientation)

        eta  =   s.     field(Coord(0,0,0)).E1().real   \
               / s_open.field(Coord(0,0,0)).E1().real 

        eta_OK = 1.13004085163

        print(eta, "expected", eta_OK)

        eta_pass = abs((eta - eta_OK) / eta_OK) < eps.testing_eps

        free_tmps()

        set_circ_PML(0)
       
        self.assertTrue(eta_pass)

suite = unittest.makeSuite(SpE, 'test')        

if __name__ == "__main__":
    unittest.main()
    


