
/////////////////////////////////////////////////////////////////////////////
//
// File:     polyroot.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000327
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "polyroot.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

/////////////////////////////////////////////////////////////////////////////
//
// polyroot
//
//   Wrapper around Netlib's Jenkins-Traub algorithm to calculate roots
//   of a polynomial.
//   The coefficients are ordered by decreasing powers of z.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(cpoly)(Real*,Real*,int&,Real*,Real*,int*);

vector<Complex> polyroot(const vector<Complex>& coef)
{
  // Set up variables.

  Real* coef_r = new Real[coef.size()];
  Real* coef_i = new Real[coef.size()];

  for (unsigned int i=0; i<coef.size(); i++)
  {
    coef_r[i] = coef[i].real();
    coef_i[i] = coef[i].imag();
  }
  
  int N = coef.size()-1;

  Real* root_r = new Real[N];
  Real* root_i = new Real[N];

  int error;

  // Call Fortan routine.

  F77NAME(cpoly)(coef_r,coef_i,N,root_r,root_i,&error);

  if (error)
    cout << "Warning: polyroot solver did not converge." << endl;

  delete [] coef_r; delete [] coef_i;
  delete [] root_r, delete [] root_i;

  // Return results.

  vector<Complex> results;
  for (unsigned int i=0; i<N; i++)
    results.push_back(Complex(root_r[i],root_i[i]));

  return results;
}
