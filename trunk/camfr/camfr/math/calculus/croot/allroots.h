
/////////////////////////////////////////////////////////////////////////////
//
// File:     allroots.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010322
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ALLROOTS_H
#define ALLROOTS_H

#include <vector>
#include "contour.h"

/////////////////////////////////////////////////////////////////////////////
//
// allroots
//
//   Returns all the roots in the rectangle (bottom_left,top_right) of the
//   function f.
//
//   eps, mu and max_k are numeric parameters for the patterson quadrature
//   formulas.
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> allroots
  (ComplexFunction& f, const Complex& bottom_left, const Complex& top_right,
   Real eps=1e-4, Real mu=1e-4, unsigned int max_k=4);



/////////////////////////////////////////////////////////////////////////////
//
// N_roots
//
//   Returns at least N roots of the function f, starting from the rectangle
//   (bottom_left,top_right).
//
//   If needed, the search region is expanded upwards and to the right or
//   downwards and to the right..
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {ur, dr, r} ExpandDirection;

std::vector<Complex> N_roots
  (ComplexFunction& f, unsigned int N,
   const Complex& bottom_left, const Complex& top_right,
   Real eps=1e-4, Real mu=1e-4, unsigned int max_k=4, ExpandDirection dir=ur);


#endif
