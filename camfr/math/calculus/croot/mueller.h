
/////////////////////////////////////////////////////////////////////////////
//
// File:     mueller.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010328
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef MUELLER_H
#define MUELLER_H

#include <vector>
#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// mueller
//
//   Uses the mueller algorithm to find a complex root of f based on initial
//   estimates a and b.
//   prev_zeros can contain a list of known zeros of f, which will then be
//   deflated.
//
/////////////////////////////////////////////////////////////////////////////

Complex mueller(ComplexFunction& f, const Complex& a, const Complex& b,
                Real eps=1e-14, const std::vector<Complex>* prev_zeros=0,
                int maxiter=100, bool *errorptr=0, bool verbose=false);



#endif
