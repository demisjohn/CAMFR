
/////////////////////////////////////////////////////////////////////////////
//
// File:     refsection.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20021104
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "refsection.h"

/////////////////////////////////////////////////////////////////////////////
//
// RefSection::find_modes
//  
/////////////////////////////////////////////////////////////////////////////

// TODO: check PML

void RefSection::find_modes()
{
  // Check values.

  if (global.lambda == 0)
  {
    py_error("Error: wavelength not set.");
    return;
  }
  
  if (global.N == 0)
  {
    py_error("Error: number of modes not set.");
    return;
  }

  Polarisation old_pol = global.polarisation;

  const int n = int(global.N/2);
  if (2*n != global.N)
    py_print("Warning: changing N to even number.");
  global.N = n;

  //
  // TE modes.
  //

  const Complex k2 = pow(2*pi/global.lambda * m->n(), 2);
  
  Complex kx0 = (global_section.leftwall == E_wall) ? 0.0 : pi;
  Complex ky0 = (global_slab.lowerwall   == NULL)   ? 0.0 : pi;

  Complex x_offset = (global_section.rightwall == E_wall) ? 0.0 : pi/2.;
  Complex y_offset = (global_slab.upperwall    == NULL)   ? 0.0 : pi/2.;  

  for (int ix=0; ix<n; ix++)
    for (int iy=0; iy<n; iy++)
    {
      Complex kx = (x_offset - kx0 + ix*pi)/a;
      Complex ky = (y_offset - ky0 + iy*pi)/b;

      Complex kt2 = kx*kx + ky*ky;
      if (abs(kt2) > 1e-8)
      {
        Complex kz = sqrt(k2 - kt2);

        if (real(kz) < 0) 
          kz = -kz;

        if (abs(real(kz)) < 1e-12)
          if (imag(kz) > 0)
            kz = -kz;
      }
    }
  




  // TM modes.




  // Restore globals.

  global.polarisation = old_pol;
  global.N = 2*n;
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::calc_overlap_matrices
//  
/////////////////////////////////////////////////////////////////////////////

void RefSection::calc_overlap_matrices(Section2D* profile,
                                       cMatrix* O_EE, cMatrix* O_MM,
                                       cMatrix* O_EM, cMatrix* O_zz)
{
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::RefSectionMode
//  
/////////////////////////////////////////////////////////////////////////////

RefSectionMode::RefSectionMode(Polarisation pol,   const Complex& kz,
                               const Complex& kx_, const Complex& kx0_,
                               const Complex& ky_, const Complex& ky0_)
  : SectionMode(pol, kz, NULL), kx(kx_), kx0(kx0_), ky(ky_), ky0(ky0_)
{
  A = 1.0; // TODO: normalise.
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::field
//  
/////////////////////////////////////////////////////////////////////////////

Field RefSectionMode::field(const Coord& coord) const
{
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::normalise
//  
/////////////////////////////////////////////////////////////////////////////

void RefSectionMode::normalise()
{
  // Do this in creation of mode.
}
