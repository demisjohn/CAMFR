/////////////////////////////////////////////////////////////////////////////
//
// File:     fourier.cpp
// Author:   Peter.Bienstman@UGent.be
// Date:     20050208
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "fourier.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// fourier
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier(const vector<Complex>& f, const vector<Complex>& disc, 
                int M, const Complex* d, bool extend)
{
  Complex D = d ? *d : disc.back() - disc.front();
  if (extend)
    D *= 2.0;
  Complex K = 2.*pi/D;

  cVector result(2*M+1,fortranArray);

  for (int m=-M; m<=M; m++)
  {
    Complex result_m = 0.0;

    for (unsigned int k=0; k<int(disc.size()-1); k++)
    {
      Complex factor;
      if (m==0)
      { 
        factor = disc[k+1]-disc[k];
        if (extend)
          factor *= 2.0;
      }
      else
      {
        const Complex t = -I*Real(m)*K;
        factor = (exp(t*disc[k+1])-exp(t*disc[k])) / t;
        if (extend)
          factor += (exp(-t*disc[k])-exp(-t*disc[k+1])) / t;
      }

      result_m += factor * f[k];
    }
    
    result(m+M+1) = result_m / D;
  }

  if (abs(D) < 1e-12)
    result = 0.0;
  
  return result;
}

