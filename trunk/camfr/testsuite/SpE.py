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

        set_N(60)
        set_lambda(1)
        set_circ_order(1)
        set_circ_field_type(cos_type)

        # Define waveguide and wall.

        air_m = Material(1.0)
        air = Circ(air_m(10-0.5j))

        wall = E_Wall(air)

        # Define cavities.

        d = 2.75

        top = Stack(air(d/2.) + wall)
        bot = Stack(air(d/2.) + wall)

        top_open = Stack(air(d/2.) + air(d/2.))
        bot_open = Stack(air(d/2.) + air(d/2.))

        source_pos  = Coord(0,0,0)
        orientation = Coord(1,0,0)

        # Calculate modification of spontaneous emission.
    
        cav = Cavity(top, bot)
        cav.set_source(source_pos, orientation)

        cav_open = Cavity(top_open, bot_open)
        cav_open.set_source(source_pos, orientation)

        eta  =   top.     field(Coord(0,0,0)).E1().real   \
               / top_open.field(Coord(0,0,0)).E1().real 

        eta_OK = 1.13004085163

        print eta, "expected", eta_OK

        eta_pass = abs(eta - eta_OK) < eps.testing_eps
       
        self.failUnless(eta_pass)

suite = unittest.makeSuite(SpE, 'test')        

if __name__ == "__main__":
    unittest.main()
    


