#! /usr/bin/env python

#############################################################################
#
# Surface Plasmon Biosensor
#
#############################################################################

from camfr_work import *
from cmath import *

import unittest, eps

class plasmon_biosensor(unittest.TestCase):
    def test_plasmon_biosensor(self):

        """Plasmon Biosensor"""

       
        print "Running plasmon Biosensor..."

        # Initialisation calculation parameters.
        
        set_lambda(1.550)
	set_polarisation(TM)
	set_N(20)
	set_keep_all_1D_estimates(False)
	set_solver(stretched_ASR)	
	set_mode_surplus(20)
	set_low_index_core(True)
        
        # Secondary calculation parameters.	
	  
	set_upper_wall(slab_H_wall)
	set_lower_wall(slab_H_wall)
	set_upper_PML(-0.05)
	set_lower_PML(-0.05)

	set_orthogonal(False) 
    
	# Thicknesses and Lengths
        
    	d_wav  = 0.220
	d_gold = 0.060
	length = 10.0

	# Initializing and Calculating the Slab.
        
    	buffer  = Material(1.33)
	protein = Material(1.923)
	Si      = Material(3.47640956822)
	SiO2	= Material(1.44402)
	Au      = Material(0.55653715538 - 9.93556321412*1j)

	# Definition of the Slab.
	
	sensorslab = Slab(SiO2(5.0)+Si(d_wav)+Au(d_gold)+protein(5.00-d_gold))
	inoutslab =  Slab(SiO2(5.0)+Si(d_wav)+buffer(5.0))	
  			
    	# Definition of the Stack.
        
        stack = Stack(inoutslab(0)+sensorslab(length)+inoutslab(0))
		
	# Incident Field
	inc 	= zeros(N())
	inc[0]	=  1	
	stack.set_inc_field(inc)
		
	stack.calc()		
		
        T_OK   = -25.5819930951
        T_test = 10*log10(abs(stack.T12(0,0))**2)
        
        print T_test , "expected", T_OK
        T_pass = abs((T_test-T_OK)/T_OK) < eps.testing_eps

        free_tmps()
        
	set_polarisation(TE)
	set_solver(track)
    	set_mode_surplus(1.5)
        set_low_index_core(True)

        self.failUnless(T_pass)

suite = unittest.makeSuite(plasmon_biosensor, 'test')

if __name__ == "__main__":
    unittest.main()
