#! /usr/bin/env python

###################################################################
#
# 90 deg 3dB splitter in rectangular lattice of rectangular 
# GaAs rods in air
#
###################################################################

from camfr import *

import unittest, eps

class PhC_splitter(unittest.TestCase):
    def testPhC_splitter(self):
        
        """Photonic crystal splitter"""

        print
        print "Running photonic crystal splitter..."

        set_lambda(1.5)
        set_N(50)
        set_polarisation(TE)
    
        # Set geometry parameters

        GaAs = Material(3.4)
        wg_m = Material(1.8)
        air  = Material(1.0)
  
        a = .600     # period
        r = .150/2.0 # rod radius

        set_lower_wall(slab_H_wall)

        cl = 0 # air cladding

        periods = 3  # periods above outer waveguide
        sections = 1 # intermediate 90 deg sections

        # Define slabs.

        no_rods = Slab(air(a-r+(sections+1+periods)*a+cl))

        # Central waveguide.

        cen = Slab(  air(a-r)                                       \
           + (sections+1+periods)*(GaAs(2*r) + air(a-2*r))          \
           + air(cl) )

        # Vertical section.

        ver = Slab(  air(a-r + (sections+1)*a)                      \
           + periods*(GaAs(2*r) + air(a-2*r) )                      \
           + air(cl) )

        # Outer arms.
 
        arm = Slab(  GaAs(r) + air(a-2*r)                           \
           + sections*(GaAs(2*r) + air(a-2*r))                      \
	   + air(a)                                                 \
           + periods*(GaAs(2*r) + air(a-2*r))                       \
           + air(cl) )

        # Find lowest order waveguide mode.

        wg = BlochStack(cen(2*r) + no_rods(a-2*r))
        wg.calc()

        guided = 0

        while wg.mode(guided).kz().real <= 0:
            guided += 1
        for i in arange(guided,2*N()):
            if (abs(wg.mode(i).kz().imag) < abs(wg.mode(guided).kz().imag)):
                if wg.mode(i).kz().real > 0:
                    guided = i

        guided_kz = wg.mode(guided).kz()
        guided_kz_OK = 2.88206436787+6.11678736051e-06j

        print guided_kz, "expected", guided_kz_OK

        guided_kz_pass \
          = abs((guided_kz - guided_kz_OK) / guided_kz_OK) < eps.testing_eps

        # Calculate splitter.

        splitter = Stack(  5*(cen(2*r) + no_rods(a-2*r))           \
                 +    ver(2*r) + no_rods(a-2*r)                    \
                 + 5*(arm(2*r) + no_rods(a-2*r)) )

        splitter.set_inc_field(wg.mode(guided).fw_field())
        splitter.calc()

        R = splitter.R12(0,0)
        R_OK = -0.850178703613-0.01113704677j

        print R, "expected", R_OK

        R_pass = abs((R - R_OK)/R_OK) < eps.testing_eps

        # Calculate field.

        E_field = splitter.field(Coord(a/2.,0,a/2.)).E2()
        E_field_OK = 22.450402472+55.5627710986j

        print E_field, "expected", E_field_OK

        E_field_pass = abs((E_field-E_field_OK)/E_field_OK) \
                           < 100*eps.testing_eps
        
        free_tmps()

        set_lower_wall(slab_E_wall)
        
        self.failUnless(guided_kz_pass and R_pass and E_field_pass)

suite = unittest.makeSuite(PhC_splitter, 'test')        

if __name__ == "__main__":
    unittest.main()
    
