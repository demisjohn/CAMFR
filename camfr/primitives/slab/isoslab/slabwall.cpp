
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabwall.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000927
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "../../planar/planar.h"
#include "slabwall.h"

/////////////////////////////////////////////////////////////////////////////
//
// Some useful slabs walls as global variables
//  
/////////////////////////////////////////////////////////////////////////////

const SlabWallMixed slab_E_wall (1.0,  1.0);
const SlabWallMixed slab_H_wall (1.0, -1.0);
const SlabWallMixed slab_no_wall(0.0,  1.0);



/////////////////////////////////////////////////////////////////////////////
//
// SlabWall_TBC::get_R12
//  
/////////////////////////////////////////////////////////////////////////////

Complex SlabWall_TBC::get_R12() const
{
  Complex k = 2*pi/global.lambda*m->n();
  Complex kx = sqrt(k*k - pow(Planar::get_kt(),2));

  int sign = (global.polarisation == TE) ? 1 : -1;

  return Real(sign) * (kx - kx_0) / (kx + kx_0);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabWall_TBC::get_start_field
//  
/////////////////////////////////////////////////////////////////////////////

void SlabWall_TBC::get_start_field(Complex* in, Complex* out) const
{
  Complex k = 2*pi/global.lambda*m->n();
  Complex kx = sqrt(k*k - pow(Planar::get_kt(),2));

  int sign = (global.polarisation == TE) ? 1 : -1;
  
  *in  =     kx_0 + kx;
  *out = - ( kx_0 - kx ) * Real(sign);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabWall_TBC::get_error
//  
/////////////////////////////////////////////////////////////////////////////

Complex SlabWall_TBC::get_error(const Complex& in, const Complex& out) const
{
  Complex k = 2*pi/global.lambda*m->n();
  Complex kx = sqrt(k*k - pow(Planar::get_kt(),2));

  int sign = (global.polarisation == TE) ? 1 : -1;
  
  return Real(sign)*(kx_0 - kx)*in + (kx_0 + kx)*out;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabWall_PC::get_R12
//  
/////////////////////////////////////////////////////////////////////////////

Complex SlabWall_PC::get_R12() const
{ 
  s.calcRT();

  // Build transfer matrix.

  const Complex& R12 = s.get_R12(); const Complex& R21 = s.get_R21();
  const Complex& T12 = s.get_T12(); const Complex& T21 = s.get_T21();

  cMatrix T(2,2,fortranArray);
  T(1,1) = T12 - R21*R12/T21; T(1,2) = R21/T21;
  T(2,1) =          -R12/T21; T(2,2) = 1.0/T21;

  if (abs(T21) < 1e-10)
    T = 0;
  
  // Solve eigenproblem.

  cVector alpha(2,fortranArray);
  cMatrix FB(2,2,fortranArray);

  alpha = eigenvalues(T, &FB);

  Complex R1 = FB(2,1)/FB(1,1);
  Complex R2 = FB(2,2)/FB(1,2);
  
  if (    ( abs(abs(alpha(1))-1) < 1e-3 )
       && ( abs(abs(alpha(2))-1) < 1e-3 ) )
  {
    //cout << "Warning: both eigenvalues on the unit circle. " << endl;
    return ( abs(R1) > abs(R2) ) ? R2 : R1;
  }

  return ( abs(alpha(1)) > abs(alpha(2)) ) ? R2 : R1;
}
