#! /usr/bin/env python

###################################################################
#
# Test metallic structures.
#
###################################################################

from cmath import *
from camfr import *
from Numeric import *

import unittest, eps

class metal_splitter(unittest.TestCase):
    def testRT(self):
        
        """Metal splitter"""

        print
        print "Running metal splitter..."

        set_lambda(1.5)
        set_N(60)
        set_unstable_exp_threshold(1e-8)

        # Set geometry parameters

        GaAs = Material(3.5)
        met  = Material(-10j)
        air  = Material(1.0)
  
        a = .600     # period
        r = .150/2.0 # rod radius

        set_left_wall(slab_H_wall)

        PML = 0.0

        cl = 0 # air cladding

        periods = 3  # periods above outer waveguide
        sections = 1 # intermediate 90 deg sections

        # Define slabs.

        no_rods = Slab(air(a-r+(sections+1+periods)*a+cl-PML*1j))

        inc_wg = Slab(GaAs(a/2.)+ air(no_rods.width()-a/2.))

        # Central waveguide.
 
        cen = Slab(  air(a-r)                                        \
           + (sections+1+periods)*(met(2*r) + air(a-2*r))            \
           + air(cl-PML*1j) )

        # Outer arms.
 
        arm = Slab(  met(r) + air(a-2*r)                             \
           + sections*(met(2*r) + air(a-2*r))                        \
	   + air(a)                                                  \
           + periods*(met(2*r) + air(a-2*r))                         \
           + air(cl-PML*1j) )

        # Vertical section.

        ver = Slab(  air(a-r + (sections+1)*a)                       \
           + periods*(met(2*r) + air(a-2*r) )                        \
           + air(cl-PML*1j) )

        # Calculate splitter.

        splitter = Stack(inc_wg(3*a) + 5*(cen(2*r) + no_rods(a-2*r)) \
                 +    ver(2*r) + no_rods(a-2*r)                      \
                 + 5*(arm(2*r) + no_rods(a-2*r)) )


        splitter.calc()

        R = splitter.R12(0,0)
        R_OK = 0.784056622162+0.388717728367j

        print R, "expected", R_OK
        
        R_pass = abs((R - R_OK)/R_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(R_pass)

suite = unittest.makeSuite(metal_splitter, 'test')        

if __name__ == "__main__":
    unittest.main()
