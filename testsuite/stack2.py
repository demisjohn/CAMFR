#! /usr/bin/env python

####################################################################
#
# stack with zero length first slab
#
####################################################################

from camfr_work import *

import unittest, eps

class stack2(unittest.TestCase):
    def teststack2(self):
        
        """Stack 2"""

        print
        print "Running stack 2..."

        set_N(10)
        set_lambda(1/.2)

        # Define materials.
        
        m = Material(3.37)
        a = Material(1.0)
        
        # Define geometry.
        
        set_lower_wall(slab_H_wall)
        set_upper_PML(-0.05)

        slab1 = Slab(a(1.26419)+m(0.271616)+a(0.728384)+m(0.271616)+\
                a(0.728384)+m(0.271616)+a(0.728384)+m(0.271616)+a(1.46419))
    
        slab2 = Slab(m(0.0607762)+a(1.16261)+m(0.353235)+a(0.646765)+\
                m(0.353235)+a(0.646765)+m(0.353235)+a(0.646765)+m(0.353235)+\
                a(1.42338-5.19999999948e-06))

        s = Stack(slab1(0)+slab2(0))
        s.calc()

        R = s.R12(0,0)
        R_OK = -0.0579512119293+0.00108957290649j
        
        print R, "expected", R_OK
        R_pass = abs((R - R_OK) / R_OK) < eps.testing_eps

        free_tmps()

        set_lower_wall(slab_E_wall)
        set_upper_PML(0)
        
        self.failUnless(R_pass)

suite = unittest.makeSuite(stack2, 'test')        

if __name__ == "__main__":
    unittest.main()
