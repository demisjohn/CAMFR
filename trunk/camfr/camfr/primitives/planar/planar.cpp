
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

  mode = new PlanarMode(TE, 0.0, 0.0, 0.0, 0.0, *this);
  calc_kz();
}



/////////////////////////////////////////////////////////////////////////////
//
// Planar::set_theta
//
/////////////////////////////////////////////////////////////////////////////

void Planar::set_theta(Real theta_radians)
{  
  if (global.lambda == 0.0)
    py_error("Error: wavelength not set.");

  kt = (2*pi / global.lambda * core->n()) * sin(theta_radians);
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



/////////////////////////////////////////////////////////////////////////////
//
// PlanarMode::field
//
//   Note that the z-dependence is taken care of elsewhere.
//
/////////////////////////////////////////////////////////////////////////////

Field PlanarMode::field(const Coord& coord) const
{
  Complex k0 = 2*pi/global.lambda;
  Complex kt = geom->get_kt();
  geom->calc_kz();

  Field field;

  if (pol == TE)
  {
    const Complex C = 1.0 / (k0*c) * B;
        
    field.E1 = 0.0;
    field.E2 = 1.0;
    field.Ez = 0.0;
    
    field.H1 = -C * kt * exp(-I*kt*coord.c1);
    field.H2 = 0.0;
    field.Hz =  C * kz;
  }
  else
  {
    const Complex C = 1.0 / (k0*c) * A;
 
    field.H1 = 0.0;
    field.H2 = 1.0;
    field.Hz = 0.0;

    field.E1 =  C * kt * exp(-I*kt*coord.c1);
    field.E2 = 0.0;
    field.Ez = -C * kz;
  }
  
  return field;
}
