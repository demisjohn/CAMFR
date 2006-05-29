
/////////////////////////////////////////////////////////////////////////////
//
// File:     material.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990407
// Version:  1.01
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#include "material.h"

/////////////////////////////////////////////////////////////////////////////
//
// BaseMaterial::operator()
//
/////////////////////////////////////////////////////////////////////////////

const Material_length BaseMaterial::operator()(const Complex& d) const
{
  if (real(d) < 0)
    py_error("Warning: negative real length of material.");
  
  return Material_length(const_cast<BaseMaterial*>(this),d);
}



/////////////////////////////////////////////////////////////////////////////
//
// Material::set_epsr_mur
//
/////////////////////////////////////////////////////////////////////////////

void Material::set_epsr_mur(const Complex& epsr_, const Complex& mur_)
{
  if ((real(epsr_) > 0.0) && (real(mur_) > 0.0)) // Right-handed materials.
  {
    i_n    = sqrt(epsr_*mur_);
    i_etar = sqrt(epsr_/mur_);
  }

  else // Use safe square roots for left-handed materials.
  {
    i_n    = sqrt(abs(epsr_)*abs(mur_))*exp(I*(arg(epsr_)+arg(mur_))/2.0);
    i_etar = sqrt(abs(epsr_)/abs(mur_))*exp(I*(arg(epsr_)-arg(mur_))/2.0);
  };
}




// Note: following is obsolete and will be removed

/////////////////////////////////////////////////////////////////////////////
//
// signedsqrt(kz2,n)
//   take sqrt of kz^2, but make sure the propagation constant
//   has the right sign (exp(jwt) time convention)
//
/////////////////////////////////////////////////////////////////////////////

/*
Complex signedsqrt(const Complex& kz2, const Material& material)
{
  const Complex n = material.n();
  
  Complex kz = sqrt(kz2);    
  if (real(kz) < 0.0)
    kz = -kz; // forward propagating wave
  
  return kz;
  
  // the previous check suffices when real(kz) != 0

  if (abs(imag(kz)) < 1e-6)
    return kz;

  int sign = 1;
  if (imag(n) > 0.0) // gain, so kz.im should be pos.
  {  
    if (imag(kz) < 0.0)
      sign = -1;
  }
  else               // loss, so kz.im should be neg.
  {
    if (imag(kz) > 0.0)
      sign = -1;
  }
  
  return sign*kz;
  
} 

*/
