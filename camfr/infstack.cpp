
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

#include "infstack.h"

// TODO: diag scatterers

/////////////////////////////////////////////////////////////////////////////
//
// InfStack::InfStack
//
/////////////////////////////////////////////////////////////////////////////

InfStack::InfStack(const Expression& e, const Complex& W_=0)
  : s(e), W(W_) 
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
    const BlochMode* mode = dynamic_cast<const BlochMode*>(s.get_mode(k));

    if (abs(mode->get_kz()) < 1e-5) // Detect presence of dummy mode.
    {
      j++;
      continue;
    }
    
    // Determine if the mode is forward or backward.

    bool fw = false;
    
    if (imag(mode->get_kz()) < -1e-3)
      fw = true;
    else if (imag(mode->get_kz()) > 1e-3)
      fw = false;
    else
    {
      fw = real(mode->get_kz()) > 0.0;

      if (fw && (mode->S_flux(0,real(W), 1e-5) < 0))
        cout << "Unexpected power flow for mode " << mode->get_kz() << endl;
    }

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

  // Calculate R12. Set T12 to the unity matrix for practical purposes.

  allocRT();

  R12.reference(multiply(B, invert_svd(F)));
  // TODO: switch to LQ factorisation.

  T12 = 0.0;
  for (int i=1; i<=N; i++)
    T12(i,i) = 1.0;

  R21 = 0.0;
  T21 = 0.0;
}
