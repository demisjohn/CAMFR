
/////////////////////////////////////////////////////////////////////////////
//
// File:     fourier.h
// Author:   Peter.Bienstman@UGent.be
// Date:     20050208
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include "../../linalg/linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// fourier
//
//  Returns fourier expansion of order m = -M,...,M of a piecewise 
//  continuous function f given by:
//
//    disc[0] -> disc[1]   : f[0]
//    ...
//    disc[n] -> disc[n+1] : f[n]
//
//  Basis functions are exp(j.m.2.pi/d.x). (d can be overridden).
//
//  Also has the option to extend the profile by adding mirror image at x=0.
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier(const std::vector<Complex>& f, 
                const std::vector<Complex>& disc, int M,
                const Complex* d=0, bool extend=false);


#endif
