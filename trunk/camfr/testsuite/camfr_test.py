#!/usr/bin/env python

##############################################################################
#
# Regression tests for CAMFR 
#
##############################################################################

import unittest

import blazed_grating, substacks, planarTE, planarTM, VCSEL, SpE, \
       PhC_splitter, expressions, cladding, grating, coupled, \
       degenerate, precision, infstack, taper, rods, \
       TEM_field, grating2, polariton, gaussian, wg, dent, stack0, stack1

alltests = unittest.TestSuite((blazed_grating.suite, substacks.suite, 
       planarTE.suite, planarTM.suite, VCSEL.suite, SpE.suite,
       PhC_splitter.suite, expressions.suite, coupled.suite,
       cladding.suite, grating.suite, degenerate.suite, stack0.suite,
       precision.suite, infstack.suite, taper.suite, rods.suite, wg.suite,
       TEM_field.suite, grating2.suite, polariton.suite, gaussian.suite,
       dent.suite, stack1.suite))

if __name__ == "__main__":
    r = unittest.TextTestRunner()
    r.run(alltests)
