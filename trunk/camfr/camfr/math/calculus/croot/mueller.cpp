
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

#include <iostream>
using namespace std;

#include <sstream>
#include "mueller.h"
#include "../../../util/vectorutil.h"

#ifdef _WIN32
#include <float.h>
#define ISNAN _isnan
#else
#define ISNAN isnan
#endif



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
 
    if ( (abs(corr) > 1e6*abs(z2)) || ISNAN(abs(z1)) )
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

  if ( (n == maxiter) || ISNAN(abs(z1)) )
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

    // Less general, assumes f is even.
    
    if (abs(new_root+z0) < eps)
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
      Complex new_root = mueller(f, z0[i]+0.001, z0[i]+.001*I+0.001, eps,
                                 &deflate, maxiter_i, &error);

      if (abs(new_root-z0[i]) > 0.1)
        error = true;

      // Less general, assumes f is even.

      if (abs(new_root+z0[i]) < eps)
        error = true;

      //std::cout << "candidate " << z0[i] << " ->" << new_root;  
      //std::cout << (error ? " reject " : " accept ") << std::endl;
      //std::cout << "----------" << std::endl;

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



/////////////////////////////////////////////////////////////////////////////
//
// mueller
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> mueller
  (ComplexFunction& f,const std::vector<Complex>& z0,Real eps,int maxiter,
   ComplexFunction* transform, int verbosity)
{
  // Calculate roots.

  vector<Complex> z1;
  for (unsigned int i=0; i<z0.size(); i++)
  {
    bool error = false;
    bool verbose = verbosity == 2;
    vector<Complex> deflate;
    Complex new_root = mueller(f, z0[i]+0.001, z0[i]+.001*I, eps,
                               &deflate, maxiter, &error, verbose);

    if (error)
    {
      std::ostringstream s;
      s << "Mueller solver failed to converge for ";
      if (!transform)
        s << z0[i] << ".";
      else
        s << (*transform)(z0[i]) << ".";      
      py_print(s.str());
    }
    else
      z1.push_back(new_root);

    if (verbosity > 0)
    {
      std::ostringstream s;
      if (!transform)
        s << i << " " << z0[i] << " --> " << new_root;
      else
        s << i << " " << (*transform)(z0[i]) << " --> "
          << (*transform)(new_root);
    
      py_print(s.str());
    }
  }

  // Dectect doubles.

  vector<Complex> z_final;
  vector<vector<Complex> > duplicates;
  
  vector<bool> dealt_with;
  for (unsigned int i=0; i<z1.size(); i++)
    dealt_with.push_back(false);
  
  for (unsigned int i=0; i<z1.size(); i++)
  {
    if (dealt_with[i])
      continue;
   
    int occurrence = 1;
    for (int j=i+1; j<z1.size(); j++)
      if (abs(z1[j] - z1[i]) < 1e-6)
      {
        occurrence++;
        dealt_with[j] = true;
      }
    
    if (occurrence == 1)
      z_final.push_back(z1[i]);
    else
    {
      vector<Complex> duplicates_i;
      for (int k=0; k<occurrence; k++)
        duplicates_i.push_back(z1[i]);
      
      duplicates.push_back(duplicates_i);
    }
  }

  // Deflate doubles.

  for (unsigned int i=0; i<duplicates.size(); i++)
  {
    const Complex z_cluster = duplicates[i][0];
    
    // Deflate neighbouring zeros.
    
    vector<Complex> deflate;
 
    for (unsigned int j=0; j<z_final.size(); j++)
      if (    (abs(z_final[j]-z_cluster) < 0.05)
           && (abs(z_final[j]-z_cluster) > 0.001) )
      {
        //std::cout << "Deflating other "<<i<< " " <<(*transform)(z_final[j]) 
        //          << z_final[j] << std::endl;
        //deflate.push_back(z_final[j]);
      }

    z_final.push_back(z_cluster);
    
    for (unsigned int j=1; j<duplicates[i].size(); j++)
    {
      deflate.push_back(z_cluster);
      std::cout << "Deflating cluster " << i << " " << (*transform)(z_cluster)
                << z_cluster << std::endl;
      
      bool error = false;
      bool verbose = verbosity == 2; 
      int maxiter_i = maxiter * deflate.size();
      Complex new_root = mueller(f, z_cluster+0.001, z_cluster+.001*I, eps,
                                 &deflate, maxiter_i, &error, verbose);

      if (error)
        py_print("Mueller solver failed to converge.");
      else if (is_present(z1, new_root, 1e-10))
        py_print("Discarding possibly bad estimate.");
      else
        z_final.push_back(new_root);

      if (verbosity > 0)
      {
        std::ostringstream s;
        if (!transform)
          s << i << " " << z0[i] << " --> " << new_root;
        else
          s << i << " " << (*transform)(z_cluster) << " --> "
            << (*transform)(new_root);
    
        py_print(s.str());
      }
    }
  }

  return z_final;
}

