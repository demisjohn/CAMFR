
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson.cpp
// Author:   Peter.Bienstman@rug.ac.be
//           converted to C++ from the CACM algorithms
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "patterson.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// patterson
//
//   This quadrature program uses formulae due to T. N. L. Patterson,
//   Mathematics of computation, Volume 22, 1968, pages 847-856, as
//   modified by F. T. Krogh and W. V. Snyder, ACM Transactions on
//   Mathematical Software 17, 4 (December 1991) pp 457-461.  It is a
//   functional replacement for Algorithm 468, T. N. L. Patterson,
//   Communications of the ACM 16, 11 (November 1973) 694-699.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson(RealFunction& f, Real a, Real b, Real eps,
               bool* error_ptr, unsigned int max_k,
               Real* abs_error)
{
  // Check if a and b are different.

  if (1. + abs(a-b) <= 1.)
    return 0.0;
  
  // Check and coerce k.
  
  if (max_k < 2)
  {
    py_print("Warning: increasing max_k to 2.");
    max_k = 3;
  }
    
  if (max_k > 8)
  {
    py_print("Warning: restricting max_k to 8.");
    max_k = 8;
  }
  
  // Include the array 'p' with the coefficients used in these formulas. 

  #include "patterson_coeff.cpp"

  // Define constants and workspace containing previous function evaluations.

  const Real diff = 0.5*(b-a);
  
  Real work[18];
  
  // Apply 1-point Gauss formula (midpoint rule).
  
  Real fx = f(a+diff); // Don't write 0.5*(b+a) if radix of arithmetic != 2.
  
  work[1] = fx;
  Real acum = fx*(b-a);
  Real result = acum;
  
  // Go on to the next formulas.

  // fl, fh contain subscripts to indicate where to store integrand
  //        samples in 'work'.
  // kl, kh are bounds for indexes at which integrand samples are
  //        retrieved from 'work' in order to begin applying a quadrature
  //        formula.
  // kx     is a list of bounds of subscripts into kh and kl.  kh is
  //        indexed by k.  The set of indices from which to retrieve
  //        integrand samples is given by the set of bounds in kl and kh
  //        indexed by kx[k-1]+1 to kx[k] inclusively.

  static const int fl[] = {0, 0, 2, 3, 5, 9,12,14, 1};
  static const int fh[] = {0, 0, 2, 4, 8,16,17,17, 0};
  static const int kl[] = {0,       1, 1, 1, 1, 1, 3, 5, 9, 5, 9,12};
  static const int kh[] = {0,       1, 2, 4, 8,16, 3, 6,17, 5, 9,17};
  static const int kx[] = {0, 0,    1, 2, 3, 4, 5,       8,      11};

  int ip=1; // Index in coefficient array 'p'.
  int jh=0;

  Real prev_result;
  
  for (int k=2; k<=max_k; k++)
  {
    prev_result = result;
    Real prev_acum  = acum;
    
    acum = 0.0;
    
    // Compute contribution to current estimate due to function
    // values used in previous formulas.
    
    for (int kk=kx[k-1]+1; kk<=kx[k]; kk++)
      for (int j=kl[kk]; j<=kh[kk]; j++)
        acum += p[ip++] * work[j];

    // Compute contribution from new function values.

    int jl=jh+1;  
        jh=jl+jl-1;
    int j1=fl[k];
    int j2=fh[k];

    for (int j=jl; j<=jh; j++)
    {
      Real x = p[ip++]*diff;   
      fx = f(a+x)+f(b-x);
      
      acum += p[ip++]*fx;

      if (j1 <= j2)
        work[j1++] = fx;
    }

    acum = diff*acum + 0.5*prev_acum;
    result = acum;
    
    if (abs(result-prev_result) <= abs(eps*result))
    {
      if (error_ptr)
        *error_ptr = false;

      if (abs_error)
        *abs_error = result - prev_result;
      
      return result;
    }
  }

  // Failed to converge.

  if (error_ptr)
    *error_ptr = true;

  if (abs_error)
    *abs_error = result - prev_result;

  return result;
}
