
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
//  Create two other estimates by varying real and imaginary part of
//  initial estimate until the function values are different enough.
//
/////////////////////////////////////////////////////////////////////////////

Complex mueller(ComplexFunction& f, const Complex& a,
                Real eps, const vector<Complex>* prev_zeros,
                int maxiter, bool *errorptr, bool verbose)
{
  const Real min_difference = 1e-6;
  
  const Complex fa = f(a);
  
  Complex b = a;
  Real increment = eps;
  do
  {
    b += increment;
    increment *= 2.0;
  }
  while (abs((f(b) - fa)/fa) < min_difference);


  Complex c = a;
  increment = eps ;
  do
  {
    c += increment*I;
    increment *= 2.0;
  }
  while (abs((f(c) - fa)/fa) < min_difference);
  
  return mueller(f, a, b, c, eps, prev_zeros, maxiter, errorptr, verbose);
}



/////////////////////////////////////////////////////////////////////////////
//
// mueller
//
//  Supply two estimates, and create the third from Newton's method.
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
      py_print("The two initial estimates are too close together.");

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

  Complex z1 = z2 - f2 * (z2 - z3) / (f2 - f3);

  // Detect convergence problems.

  if (1. + abs(f2-f3) <= 1.)
  {
    if (verbose)
      cout << "The two initial estimates are too close together." << endl;

    //if (errorptr)
    //  *errorptr = true;
    
    return a;
  }
  
  if (abs(z1-z2) > 1e4*abs(z2))
  {
    if (verbose)
      cout << "Convergence problems. Possible bad initial estimates." << endl;

    if (errorptr)
      *errorptr = true;
    
    return a;
  }

  // Call classic Mueller solver.

  return mueller(f, a, b, z1, eps, prev_zeros, maxiter, errorptr, verbose);
  
}



/////////////////////////////////////////////////////////////////////////////
//
// mueller
//
//  'Classic' version, with three initial estimates.
//
/////////////////////////////////////////////////////////////////////////////

Complex mueller(ComplexFunction& f, const Complex& a, const Complex& b,
                const Complex& c,
                Real eps, const vector<Complex>* prev_zeros,
                int maxiter, bool *errorptr, bool verbose)
{
  Complex z3 =   a;  Complex z2 =   b;  Complex z1 = c;
  Complex f3 = f(a); Complex f2 = f(b);

  if (prev_zeros) // Deflate them.
    for (unsigned int i=0; i<prev_zeros->size(); i++)
    {
      f3 /= (a - (*prev_zeros)[i]);
      f2 /= (b - (*prev_zeros)[i]);      
    }

  Complex dz23 = (f2 - f3) / (z2 - z3);

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
      if (verbose)
        cout << "Error increased dramatically." << endl;
      
      converged = false;
      break;
    }
    
    // Avoid division by zero in case convergence has happened.

    if ( (abs(z1-z2) < 1e-15) || (abs(z1-z3) < 1e-15) )
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
      cout << "Returned estimate  = " << z1 << endl;
    }

  if (errorptr)
    *errorptr = !converged;  
  
  return z1;
}



/////////////////////////////////////////////////////////////////////////////
//
// mueller_multiple
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> mueller_multiple
  (ComplexFunction& f, const Complex& z0, Real eps,
   const std::vector<Complex>* prev_zeros, int maxiter)
{
  vector<Complex> roots;
  roots.push_back(z0);

  vector<Complex> deflate;
  deflate.push_back(z0);
  if (prev_zeros)
    deflate.insert(deflate.end(), prev_zeros->begin(), prev_zeros->end());

  bool error = false;
  do
  {
    Complex new_root = mueller(f, z0+0.001, z0+.001*I, eps, 
                               &deflate, maxiter, &error);

    if (abs(new_root-z0) > 0.1)
      error = true;

    if (!error)
    {
      roots.push_back(new_root);
      deflate.push_back(new_root);
    }
    
  }
  while (!error);
  
  return roots;
}



/////////////////////////////////////////////////////////////////////////////
//
// mueller_multiple
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> mueller_multiple
 (ComplexFunction& f, const std::vector<Complex>& z0, Real eps,
  const std::vector<Complex>* prev_zeros, int maxiter)
{
  // Create clustered vector of zeros.

  vector<vector<Complex> > roots;
  for (unsigned int i=0; i<z0.size(); i++)
  {
    vector<Complex> roots_i;
    roots_i.push_back(z0[i]);
    roots.push_back(roots_i);
  }
  
  for (unsigned int i=0; i<z0.size(); i++)
  {
    // Determine which zeros to deflate.

    vector<Complex> deflate;

    if (prev_zeros)
      deflate.insert(deflate.end(), prev_zeros->begin(), prev_zeros->end());
 
    for (unsigned int j=0; j<z0.size(); j++)
      if (abs(z0[j]-z0[i]) < 0.5)
        deflate.insert(deflate.end(), roots[j].begin(), roots[j].end());    

    // Look for degenerate zeros.
    
    int maxiter_i = maxiter * deflate.size();
    bool error = false;
    do
    {      
      Complex new_root = mueller(f, z0[i]+0.001, z0[i]+.001*I, eps,
                                 &deflate, maxiter_i, &error);

      if (abs(new_root-z0[i]) > 0.1)
        error = true;

      if (!error)
      {
        deflate.push_back(new_root);
        roots[i].push_back(new_root);
        maxiter_i += maxiter_i;
      }
    }
    while (!error);
  }

  // Return all zeros.

  vector<Complex> allroots;

  for (unsigned int i=0; i<roots.size(); i++)
    allroots.insert(allroots.end(), roots[i].begin(), roots[i].end());

  return allroots;
}


