
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
//   Uses the mueller algorithm to find a complex root of f based on 
//   1, 2 or 3 initial estimates.
//   prev_zeros can contain a list of known zeros of f, which will then be
//   deflated.
//
/////////////////////////////////////////////////////////////////////////////

Complex mueller(ComplexFunction& f, const Complex& a,
                Real eps=1e-14, const std::vector<Complex>* prev_zeros=0,
                int maxiter=100, bool *errorptr=0, bool verbose=false);

Complex mueller(ComplexFunction& f, const Complex& a, const Complex& b,
                Real eps=1e-14, const std::vector<Complex>* prev_zeros=0,
                int maxiter=100, bool *errorptr=0, bool verbose=false);

Complex mueller(ComplexFunction& f, const Complex& a, const Complex& b,
                const Complex& c,
                Real eps=1e-14, const std::vector<Complex>* prev_zeros=0,
                int maxiter=100, bool *errorptr=0, bool verbose=false);



/////////////////////////////////////////////////////////////////////////////
//
// mueller_multiple
//
//  Given an inital zero z0, locate additional zeros at or around z0.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> mueller_multiple
  (ComplexFunction& f, const Complex& z0, Real eps=1e-14, 
   const std::vector<Complex>* prev_zeros=0, int maxiter=50);

std::vector<Complex> mueller_multiple
  (ComplexFunction& f, const std::vector<Complex>& z0, Real eps=1e-14, 
   const std::vector<Complex>* prev_zeros=0, int maxiter=50);




#endif
