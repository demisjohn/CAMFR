
/////////////////////////////////////////////////////////////////////////////
//
// File:     planar.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000111
// Version:  1.2
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "planar.h"

/////////////////////////////////////////////////////////////////////////////
//
// Planar::kt
//
/////////////////////////////////////////////////////////////////////////////

Complex Planar::kt = 0.0;



/////////////////////////////////////////////////////////////////////////////
//
// Planar::Planar
//
/////////////////////////////////////////////////////////////////////////////

Planar::Planar(Material& m) : MonoWaveguide(&m)
{ 
  // Create mode and set E_cst and H_cst to get correct results in the
  // generalised Fresnel formulas (interface.cpp).
  // 'mode' is deleted in ~MonoWaveguide.

  Complex kz = calc_kz();
  
  mode = (global.polarisation == TE)
    ? new Mode(TE, kz, -kz,        0.0,        1.0/core->mur() )
    : new Mode(TM, kz, -kz, 1.0/core->epsr(),     0.0);
}



/////////////////////////////////////////////////////////////////////////////
//
// Planar::calc_kz()
//
/////////////////////////////////////////////////////////////////////////////

Complex Planar::calc_kz()
{
  Complex k  = 2.0*pi / global.lambda * core->n() * sqrt(core->mur());
  Complex kz = sqrt(k*k - kt*kt);

  if (abs(real(kz)) > 1e-12) // If propagating, always choose forward one.
  {
    if (real(kz) < 0)
      kz = -kz;
  }
  else // On imaginary axis, choose damped one.
  {
    if (imag(kz) > 0)
      kz = -kz;
  }

  return kz;
}
