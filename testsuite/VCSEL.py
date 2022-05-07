#!/usr/bin/env python

####################################################################
#
# Finds a laser mode in a VCSEL (from COST 268 modelling exercise.)
#
####################################################################

from camfr import *

import unittest, eps

class VCSEL(unittest.TestCase):
    def testRT(self):
        
        """VCSEL"""

        print("")
        print("Running VCSEL...")

        set_lambda(.980)
        set_N(100)
        set_circ_order(1)
        set_sweep_from_previous(1)
        set_circ_PML(-0.1)

        # Define materials.

        GaAs_m   = Material(3.53)
        AlGaAs_m = Material(3.08)
        AlAs_m   = Material(2.95)
        AlOx_m   = Material(1.60)
        air_m    = Material(1.00)

        gain_m = Material(3.53)
        loss_m = Material(3.53 - 0.01j)

        set_gain_material(gain_m)

        # Define geometry parameters

        r = 4.0
        d_cladding = 4.0
 
        d_GaAs   = .06949
        d_AlGaAs = .07963

        # Define cross sections

        GaAs   = Circ(  GaAs_m(r+d_cladding))
        AlGaAs = Circ(AlGaAs_m(r+d_cladding))
        air    = Circ(   air_m(r+d_cladding))

        ox = Circ(AlAs_m(r) + AlOx_m(d_cladding))
        QW = Circ(gain_m(r) + loss_m(d_cladding))

        # Set oxide window position.

        position = 4 # 1: node, 5: antinode
        x = (5. - position) * d_AlGaAs/5.

        # Define top half of cavity.
  
        top = Stack( (GaAs(0) + AlGaAs(x)) + ox(.2*d_AlGaAs)           \
            + (AlGaAs(.8*d_AlGaAs - x) + GaAs(d_GaAs)                  \
            + 24*(AlGaAs(d_AlGaAs) + GaAs(d_GaAs)) + air(0)) )
  
        # Define bottom half of cavity.

        bottom = Stack(GaAs(.13659) + QW(.00500)                       \
          + (GaAs(.13659) + 30*(AlGaAs(d_AlGaAs) + GaAs(d_GaAs) + GaAs(0))) )
  
        # Define cavity and find laser mode.

        cavity = Cavity(bottom, top)

        cavity.find_mode(.9803, .9804)

        wavelength = get_lambda()
        gain = gain_m.gain()
        
        wavelength_OK = 0.980365834012
        gain_OK = 1182.17331249

        print(wavelength, "expected", wavelength_OK)
        print(gain, "expected", gain_OK)

        wavelength_pass \
          = abs((wavelength - wavelength_OK)/wavelength_OK) < eps.testing_eps
        gain_pass = abs((gain - gain_OK)/gain_OK) < 20*eps.testing_eps

        free_tmps()

        set_sweep_from_previous(0)

        set_circ_PML(0)

        self.assertTrue(wavelength_pass and gain_pass)

suite = unittest.makeSuite(VCSEL, 'test')        

if __name__ == "__main__":
    unittest.main()
    

