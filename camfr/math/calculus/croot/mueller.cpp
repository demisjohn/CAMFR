
/////////////////////////////////////////////////////////////////////////////
//
// File:     mueller.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010328
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "mueller.h"

#include <iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// mueller
//
/////////////////////////////////////////////////////////////////////////////

Complex mueller(ComplexFunction& f, const Complex& a, const Complex& b,
                Real eps, const vector<Complex>* prev_zeros,
                int maxiter, bool *errorptr, bool verbose)
{
  // Check input.
  
  if (1. + abs(a-b) <= 1.)
  {
    if (verbose)
      cout << "The two initial estimates are too close together." << endl;

    if (errorptr)
      *errorptr = true;
    
    return a;
  }

  // Create third estimate using Newton's method.
  
  Complex z3 =   a;  Complex z2 =   b; 
  Complex f3 = f(a); Complex f2 = f(b);

  if (prev_zeros) // Deflate them.
    for (unsigned int i=0; i<prev_zeros->size(); i++)
    {
      f3 /= (a - (*prev_zeros)[i]);
      f2 /= (b - (*prev_zeros)[i]);      
    }

  if (verbose)
  {
    cout << "Initial value of a = " << z3 << " f(a) = " << f3 << endl;
    cout << "Initial value of b = " << z2 << " f(b) = " << f2 << endl;
  }
  
  Complex dz23 = (f2 - f3) / (z2 - z3);
  Complex z1 = z2 - f2 / dz23;

  // Detect convergence problems.
  
  if ( (abs(z1-z2) > 1e4*abs(z2)) || (1. + abs(f2-f3) <= 1.) )
  {
    if (verbose)
      cout << "Convergence problems. Possible bad initial estimates." << endl;

    if (errorptr)
      *errorptr = true;
    
    return a;
  }

  // Iterate.

  bool converged = false;
  int n;
  for (n=1; n<maxiter; n++)
  {   
    Complex f1 = f(z1);

    if (prev_zeros) // Deflate them.
      for (unsigned int i=0; i<prev_zeros->size(); i++)
        f1 /= (z1 - (*prev_zeros)[i]);

    if (verbose)
      cout << "n = " << n << " z = " << z1 << " f(z) = " << f1 << endl;

    // If the error increases dramatically, we won't converge.
    
    if (abs(f1) > 1e4*(abs(f2)+abs(f3)))
    {
      converged = false;
      break;
    }
    
    // Avoid division by zero in case convergence has happened.

    if (    (1. + abs(z1-z2) <= 1.)
         || (1. + abs(z1-z3) <= 1.)
         || (1. + abs(f1-f2) <= 1.) )
    {
      converged = true;
      break;
    }
    
    // Calculate new estimate.
    
    Complex dz12  = ( f1  -  f2 ) / (z1 - z2);
    Complex dz123 = (dz12 - dz23) / (z1 - z3);
    
    Complex W = (z1 - z2)*dz123 + dz12;
    Complex w = sqrt(W*W - 4.*f1*dz123);
    Complex corr = (abs(W+w) > abs(W-w)) ? 2.*f1/(W+w) : 2.*f1/(W-w);

    z3 = z2; z2 = z1; z1 = z1 - corr;
    f3 = f2; f2 = f1;
    dz23 = dz12;

    // Check for convergence.

    if (abs(corr) > 1e6*abs(z2))
    {
      converged = false;
      break;
    }
    
    if (abs(corr) < eps*abs(z1))
    {
      converged = true;
      break;
    }
  }

  if (n == maxiter)
    converged = false;

  // Print diagnostics and return result.

  if (verbose)
    if (converged)
    {
      cout << "The relative accuracy " << eps << " has been reached." << endl;
      cout << "Number of iteration steps : " << n << endl;
      cout << "z1 = " << z1 << endl;
    }
    else
    {
      cout << "The Mueller solver failed to converge." << endl;
      cout << "Initial value of a = " << a << endl;
      cout << "Initial value of b = " << b << endl;
      cout << "Result             = " << z1 << endl;
    }

  if (errorptr)
    *errorptr = !converged;  
  
  return z1;
}
