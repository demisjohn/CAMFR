
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
  // Create dummy mode and update it in calc_kz.

  mode = new Mode(TE, 0.0, 0.0, 0.0, 0.0);
  calc_kz();
}



/////////////////////////////////////////////////////////////////////////////
//
// Planar::calc_kz()
//
/////////////////////////////////////////////////////////////////////////////

Complex Planar::calc_kz() const
{
  // Calculate kz.

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

  // Update mode and set E_cst and H_cst to get correct results in the
  // generalised Fresnel formulas (interface.cpp).

  mode->pol   = global.polarisation;
  mode->kz    =  kz; 
  mode->kz_bw = -kz;

  if (global.polarisation == TE)
  {
    mode->A = 0.0;
    mode->B = 1.0/core->mur();
  }
  else
  {
    mode->A = 1.0/core->epsr();
    mode->B = 0.0;
  }

  // Return result.

  return kz;
}
