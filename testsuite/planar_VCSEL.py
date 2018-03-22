#!/usr/bin/env python

####################################################################
#
# Planar VCSEL
#
####################################################################

from camfr import *

import unittest, eps

class planar_VCSEL(unittest.TestCase):
    def testRT(self):
        
        """planar VCSEL"""

        print('')
        print("Running planar VCSEL...")

        set_lambda(.98112)
        set_N(1)
        set_polarisation(TE)

        # Define materials.

        GaAs_m   = Material(3.53)
        AlGaAs_m = Material(3.08)
        AlAs_m   = Material(2.95)
        AlOx_m   = Material(1.60)
        air_m    = Material(1.00)

        gain_m = Material(3.53)
        loss_m = Material(3.53 - 0.01j)

        # Define geometry parameters

        d_cladding = 4.0
        d_GaAs     = .06949
        d_AlGaAs   = .07963
        d_well     = 0.005

        # Define cross sections

        GaAs   = Planar(  GaAs_m)
        AlGaAs = Planar(AlGaAs_m)
        air    = Planar(   air_m)
        ox_out = Planar(AlOx_m)
        ox_in  = Planar(AlAs_m)
        QW     = Planar(gain_m)

        GaAs.set_theta(0)

        # Set oxide window position.

        position = 1.0 # 1: node, 5: antinode
        x = (5. - position) * d_AlGaAs/5.

        # Define top half of cavity inside the aperture.
        
        top = Stack( (GaAs(0) + AlGaAs(x)) + ox_in(.2*d_AlGaAs)           \
                     + (AlGaAs(.8*d_AlGaAs - x) + GaAs(d_GaAs)            \
                     + 24*(AlGaAs(d_AlGaAs) + GaAs(d_GaAs)) + air(0)) )

        # Define the cavity.
        
        cavity = Stack(GaAs(.13649) + QW(d_well) + GaAs(.13649) )

        # Define bottom half of cavity.

        bottom = Stack(GaAs(0.1) + 30*(GaAs(d_GaAs) + AlGaAs(d_AlGaAs)) )

        vcsel = Stack(bottom + cavity + top)

        incL = zeros(N())
        incR = zeros(N())
        incR[0] = 1
        vcsel.set_inc_field(incL,incR)

        field = vcsel.field(Coord(0,0,bottom.length())).abs_E()
        field_OK = 31.197347844
        
        print(field, "expected", field_OK)
        
        field_pass = abs((field - field_OK)/field_OK) < eps.testing_eps

        free_tmps()

        self.failUnless(field_pass)

suite = unittest.makeSuite(planar_VCSEL, 'test')        

if __name__ == "__main__":
    unittest.main()
    

