
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

const Material_length BaseMaterial::operator()(const Complex& d=0.0) const
{
  if (real(d) < 0)
    cout << "Warning: negative real length of material." << endl;
  
  return Material_length(const_cast<BaseMaterial*>(this),d);
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
