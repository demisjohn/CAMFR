#!/usr/bin/env python

##############################################################################
#
# Regression tests for CAMFR 
#
##############################################################################

import unittest

import blazed_grating, substacks, planarTE, planarTM, VCSEL, SpE, \
       PhC_splitter, metal_splitter, expressions, cladding, grating, \
       degenerate, ADR_solver, precision, infstack, taper, rods, \
       TEM_field, grating2, polariton, gaussian

alltests = unittest.TestSuite((blazed_grating.suite, substacks.suite, 
       planarTE.suite, planarTM.suite, VCSEL.suite, SpE.suite,
       PhC_splitter.suite, metal_splitter.suite, expressions.suite,
       cladding.suite, grating.suite, degenerate.suite, ADR_solver.suite,
       precision.suite, infstack.suite, taper.suite, rods.suite,
       TEM_field.suite, grating2.suite, polariton.suite, gaussian.suite))

if __name__ == "__main__":
    r = unittest.TextTestRunner()
    r.run(alltests)
