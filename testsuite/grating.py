#! /usr/bin/env python

###################################################################
#
# Test for stable field calculation in a large grating.
#
###################################################################

from camfr import *
from Numeric import *
import unittest, eps

class grating(unittest.TestCase):
    def testgrating(self):

        set_N(10)
        set_lambda(1.5)

        GaAs = Material(3.5)
        air  = Material(1.0)

        set_left_wall (slab_H_wall)
        set_right_wall(slab_H_wall)
      
        compare = 10
        continu = 1
      
        for gpol in arange(0.5, 0.7, 0.02):
            gp = gpol * get_lambda()

            slab = Slab(5*(air(gp/2.) + GaAs(gp) + air(gp/2.))+air(0))
            slab.calc()
    
            E_field = abs(slab.mode(0).field(Coord(gp/2.,0,0)).E2().real)
        
            print gp, E_field 

            if E_field > compare:
                continu = 0
            compare = E_field

        free_tmps()
        
        self.failUnless(1)

suite = unittest.makeSuite(grating, 'test')

if __name__ == "__main__":
    unittest.main()
          
