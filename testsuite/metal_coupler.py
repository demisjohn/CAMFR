#! /usr/bin/env python

#############################################################################
#
# TM metal grating coupler
#
#############################################################################

from camfr import *
from cmath import *

import unittest, eps

class metal_coupler(unittest.TestCase):
    def test_metal_coupler(self):

        """TM metal grating coupler"""

        print('')
        print("Running TM metal grating coupler...")

        # Define calculation parameters

        set_lambda(1.55)
        set_N(100)
        set_solver(stretched_ASR)
        set_mode_surplus(5)
        set_polarisation(TM)
        set_keep_all_1D_estimates(False)
        set_low_index_core(True)
        set_degenerate(0)
        set_chunk_tracing(0)
        set_orthogonal(False)

        # Boundary conditions.

        pml = 0.5
        set_lower_PML(-pml)
        set_upper_PML(-pml)

        # Define some parameters that will be used in the calculations.

        nom             = 100
        no_of_grooves   = 12
        d_sub           = 6
        period          = 0.96
        ff              = 0.1
        metal_height    = 0.075
        d_air           = 10
        guide_thickness = 0.220
        dclad           = 2.0

        # Define materials.

        SiO2  = Material(1.447)
        Si    = Material(3.476)
        Air   = Material(1.0)
        H20   = Material(1.33)
        nfib  = 1.001
        iml   = Material(nfib)
        metal = Material(0.58-10.75*1j)

        # Define slabs.

        waveguide   = Slab(Si(d_sub) + SiO2(dclad) + Si(guide_thickness) + \
                           iml(d_air))
        metal_tooth = Slab(Si(d_sub) + SiO2(dclad) + Si(guide_thickness) + \
                           metal(metal_height) + iml(-metal_height + d_air))

        # Define stack.

        stack = Stack(waveguide(1.0) + no_of_grooves*(metal_tooth(period*ff) + \
                                    waveguide(period*(1-ff))) + waveguide(1.0))

        # Find the guided mode.

        waveguide.calc()

        # Look for the mode with the highest field component in the core.

        guided  = 0
        campmax = 0
        for k in arange(0,nom,1):
            camp = abs(waveguide.mode(k).field(Coord(d_sub+dclad+\
                                            guide_thickness/2,0,0)).E1())
            if campmax < camp:
                campmax=camp
                guided=k

        # Set input for calculating the fields.

        inc = zeros(N())
        inc[guided] = 1
        stack.set_inc_field(inc)

        # Calculate stack.

        stack.calc()

        x_up = waveguide.width()-d_air+1.5
        up = stack.lateral_S_flux(x_up)

        # Compare vertical flux (grating coupler efficiency).

        up_OK   = 0.355387646801
        up_test = up.real

        print(up_test , "expected", up_OK)
        up_pass = abs((up_test-up_OK)/up_OK) < eps.testing_eps

        free_tmps()

        self.assertTrue(up_pass)

suite = unittest.makeSuite(metal_coupler, 'test')

if __name__ == "__main__":
    unittest.main()
