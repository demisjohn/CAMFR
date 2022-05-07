#! /usr/bin/env python

#############################################################################
#
# Surface plasmon biosensor
#
#############################################################################

from camfr import *
from cmath import *

import unittest, eps

class plasmon_biosensor(unittest.TestCase):
    def test_plasmon_biosensor(self):

        """Plasmon biosensor"""

        print('')
        print("Running plasmon biosensor...")

        # Initialisation calculation parameters.
        
        set_lambda(1.550)
        set_polarisation(TM)
        set_N(20)
        set_solver(stretched_ASR)   
        set_mode_surplus(5)
        set_low_index_core(True)
        
        # Secondary calculation parameters. 
      
        set_upper_wall(slab_H_wall)
        set_lower_wall(slab_H_wall)
        set_upper_PML(-0.05)
        set_lower_PML(-0.05)

        set_orthogonal(False) 
    
        # Thicknesses and lengths.
        
        d_wav  = 0.220
        d_gold = 0.060
        length = 10.0

        # Initializing and the slabs.
        
        buffer  = Material(1.33)
        protein = Material(1.923)
        Si      = Material(3.47640956822)
        SiO2    = Material(1.44402)
        Au      = Material(0.55653715538 - 9.93556321412*1j)
    
        sensorslab = Slab(SiO2(5.0)+Si(d_wav)+Au(d_gold)+protein(5.00-d_gold))
        inoutslab =  Slab(SiO2(5.0)+Si(d_wav)+buffer(5.0))  
            
        # Definition of the stack.
        
        stack = Stack(inoutslab(0)+sensorslab(length)+inoutslab(0))
        stack.calc()
        
        T    = 10*log10(abs(stack.T12(0,0))**2) 
        T_OK = -25.5819930951
        print(T, "expected", T_OK)
        T_pass = abs((T-T_OK)/T_OK) < eps.testing_eps
        
        set_polarisation(TE)
        set_solver(track)
        set_mode_surplus(1.5)
        set_low_index_core(False)
        set_upper_wall(slab_E_wall)
        set_lower_wall(slab_E_wall)
        set_upper_PML(0)
        set_lower_PML(0)
        set_orthogonal(True)

        self.assertTrue(T_pass)

suite = unittest.makeSuite(plasmon_biosensor, 'test')

if __name__ == "__main__":
    unittest.main()
