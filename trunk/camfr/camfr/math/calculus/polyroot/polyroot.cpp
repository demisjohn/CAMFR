
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

using namespace std;

#ifndef FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#ifndef FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#ifndef FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES

#define FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE

#endif
#endif
#endif

#ifdef FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#define cpoly_F  cpoly
#endif

#ifdef FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#define cpoly_F  cpoly_
#endif

#ifdef FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES
#define cpoly_F  cpoly__
#endif



/////////////////////////////////////////////////////////////////////////////
//
// polyroot
//
//   Wrapper around Netlib's Jenkins-Traub algorithm to calculate roots
//   of a polynomial.
//   The coefficients are ordered by decreasing powers of z.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void cpoly_F(Real*,Real*,int&,Real*,Real*,int*);

vector<Complex> polyroot(const vector<Complex>& coef)
{
  // Set up variables.
  
  Real coef_r[coef.size()], coef_i[coef.size()];

  for (unsigned int i=0; i<coef.size(); i++)
  {
    coef_r[i] = coef[i].real();
    coef_i[i] = coef[i].imag();
  }
  
  int N = coef.size()-1;

  Real root_r[N], root_i[N];

  int error;

  // Call Fortan routine.

  cpoly_F(coef_r,coef_i,N,root_r,root_i,&error);

  if (error)
    cout << "Warning: polyroot solver did not converge." << endl;

  // Return results.

  vector<Complex> results;
  for (unsigned int i=0; i<N; i++)
    results.push_back(root_r[i] + I*root_i[i]);

  return results;
}
