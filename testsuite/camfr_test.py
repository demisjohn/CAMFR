#!/usr/bin/env python

##############################################################################
#
# Regression tests for CAMFR 
#
##############################################################################

import unittest

import blazed_grating, substacks, planarTE, planarTM, VCSEL, SpE, \
       PhC_splitter, metal_splitter, expressions, cladding, grating, \
       degenerate

alltests = unittest.TestSuite((blazed_grating.suite, substacks.suite, 
       planarTE.suite, planarTM.suite, VCSEL.suite, SpE.suite,
       PhC_splitter.suite, metal_splitter.suite, expressions.suite,
       cladding.suite, grating.suite, degenerate.suite))

if __name__ == "__main__":
    r = unittest.TextTestRunner()
    r.run(alltests)
