#! /usr/bin/env python

#############################################################################
#
# Surface Plasmon Polariton
#
#############################################################################

from camfr import *
from cmath import *

import unittest, eps

class surface_plasmon(unittest.TestCase):
    def test_surface_plasmon(self):

        """Surface plasmon"""

        print("")
        print("Running surface plasmon...")

        # Initialisation calculation parameters.
        
        set_polarisation(TM)
        set_lambda(1.55)
        set_solver(stretched_ASR)
        set_N(15)
        set_mode_surplus(7) 
        
        # Secondary calculation parameters. 
      
        set_low_index_core(True)
        set_keep_all_1D_estimates(False)
        
        # Gold Palik.
        
        wavelength=1.55
        n_met_r = 0.0674*wavelength*wavelength + 0.4058*wavelength - 0.2327
        n_met_im = 0.1199*wavelength*wavelength + 5.7701*wavelength + 0.5934

        # Materials.
        
        diel1 = Material(1.3)
        diel2 = Material(1.303)
        gold = Material(n_met_r -1j*(n_met_im))
    
        # Thicknesses.
        
        t_diel = 30
        t_silver = 0.02

        # Initializing and Calculating the Slab
        
        sample = Slab(diel1(t_diel)+gold(t_silver)+diel2(t_diel))
        sample.calc()   
            
        # Comparing the second mode of the structure.
        
        n_eff_OK   = 1.30342585412
        n_eff_test = sample.mode(1).n_eff().real
        
        print(n_eff_test , "expected", n_eff_OK)
        n_eff_pass = abs((n_eff_test-n_eff_OK)/n_eff_OK) < eps.testing_eps

        free_tmps()
        
        set_polarisation(TE)
        set_solver(track)
        set_mode_surplus(1.5)
        set_low_index_core(False)

        self.assertTrue(n_eff_pass)

suite = unittest.makeSuite(surface_plasmon, 'test')

if __name__ == "__main__":
    unittest.main()
