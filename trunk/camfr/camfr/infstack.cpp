
/////////////////////////////////////////////////////////////////////////////
//
// File:     infstack.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020206
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "infstack.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// InfStack::InfStack
//
/////////////////////////////////////////////////////////////////////////////

InfStack::InfStack(const Expression& e)
  : s(e)
{
  Expression e1 = 1*e; // Makes sure incidence and exit media are the same.

  inc = e1.get_inc();
  ext = e1.get_ext();
}



/////////////////////////////////////////////////////////////////////////////
//
// InfStack::calcRT()
//
/////////////////////////////////////////////////////////////////////////////

void InfStack::calcRT()
{ 
  if (!recalc_needed())
    return;
  
  // Construct F and B containing field expansion of forward modes.

  s.find_modes();

  const int N = global.N;
  
  cMatrix F(N,N,fortranArray); F = 0.0;
  cMatrix B(N,N,fortranArray); B = 0.0;

  int j=1;

  for (int k=1; k<=s.N() && j<=N; k++)
  {
    BlochMode* mode = dynamic_cast<BlochMode*>(s.get_mode(k));

    if (abs(mode->get_kz()) < 1e-5) // Detect presence of dummy mode.
    {
      j++;
      continue;
    }
    
    // Determine if the mode is forward or backward.

    bool fw = true;
    
    if (abs(imag(mode->get_kz())) < 1e-3) // Propagating mode.
    {
      Real d = mode->S_flux(0.0,real(inc->c1_size()),0.1);
      fw = (d > 0);
    }
    else // Evanescent mode.
      fw = (imag(mode->get_kz()) < 0);
        
    if (fw)
    {
      cVector f(N,fortranArray); f.reference(mode->fw_field());
      cVector b(N,fortranArray); b.reference(mode->bw_field());
      
      for (int i=1; i<=N; i++)
      {
        F(i,j) = f(i);
        B(i,j) = b(i);
      }

      j++;
    }  
  }

  // Calculate R12. Set the T's to the unity matrix for practical purposes,
  // and to avoid spurious errors.

  allocRT();

  R12.reference(multiply(B, invert_svd(F)));

  T12 = 0.0;
  for (int i=1; i<=N; i++)
    T12(i,i) = 1.0;

  R21 = 0.0;

  T21 = 0.0;
  for (int i=1; i<=N; i++)
    T21(i,i) = 1.0;

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
  last_slab_ky = global.slab_ky;
}
