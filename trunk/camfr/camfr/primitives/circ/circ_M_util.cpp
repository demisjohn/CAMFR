
/////////////////////////////////////////////////////////////////////////////
//
// File:     circ_M_util.h
// Author:   Mihai Ibanescu (michel@mit.edu)
// Date:     20020402
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "circ_M_util.h"

using std::cout;
using std::endl;

/////////////////////////////////////////////////////////////////////////////
//
// Calculates transfer matrix for interface between mat1 and mat2 at 
// a radius r.
// Uses kz as independent variable.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix transfer_matrix(const Complex& r, 
			const Material& mat1, 
			const Material& mat2,
			const Complex& kz_,
			const int& circ_order,
			const bool& scaling)
{  
  // Get indices of refraction

  const Complex n1 = mat1.n();
  const Complex n2 = mat2.n();

  if ((mat1.mur() != 1.0) || (mat2.mur() != 1.0))
  {
    cout << "Only mur=1 case implemented in transfer_matrix" << endl;
    return cMatrix();
  }

  const Complex kz = - kz_; // exp(-ikz) in CAMFR, versus exp(ikz) in Matlab

  // Set constants

  const Real    k0  = - 2*pi/global.lambda; // exp(ikz) versus exp(-ikz) 
  const Complex ord = circ_order;

  Complex k1  = sqrt(k0*k0*n1*n1 - kz*kz);
  // always put kt in the semicircle centered on the first quadrant
  if ((real(k1)+imag(k1))<0)
  {
    k1 = -k1;
    // give a warning if k1 is close to this branch-cut
    if (abs(real(k1)+imag(k1))<0.01*abs(k1))
      cout << "Warning: branch-cut in circ.cpp: " << -k1 
           << "became " << k1 << endl;
  }
  Complex k2  = sqrt(k0*k0*n2*n2 - kz*kz);
  // always put kt in the semicircle centered on the first quadrant
  if ((real(k2)+imag(k2))<0)
  {
    k2 = -k2;
    // give a warning if k2 is close to this branch-cut
    if (abs(real(k2)+imag(k2))<0.01*abs(k2))
      cout << "Warning: branch-cut in circ.cpp: " << -k2 
           << "became " << k2 << endl;
  }

  const Complex rk1 = k1*r;
  const Complex rk2 = k2*r;

  //cout << "% k1 = " << k1 << "\t k2 = " << k2 << "\t kz_ = " << kz_ << endl;
  
  const Complex c1 = (k2*n1*n1) / (k1*n2*n2);
  const Complex c2 = I * kz * ord / (k0*n2*n2*r) * (1./k2 - k2/k1/k1);
  const Complex d1 = k2/k1;
  const Complex d2 = - I * kz * ord / (k0*r) * (1./k2 - k2/k1/k1);

  // Calculate Bessel functions

  Complex H1_1, dH1_1 = dH1(circ_order, rk1, &H1_1, NULL, scaling);
  Complex H2_1, dH2_1 = dH2(circ_order, rk1, &H2_1, NULL, scaling);

  Complex H1_2, dH1_2 = dH1(circ_order, rk2, &H1_2, NULL, scaling);
  Complex H2_2, dH2_2 = dH2(circ_order, rk2, &H2_2, NULL, scaling);

  // Calculate transfer matrix

  cMatrix T(4,4,fortranArray);

  T(1,1) = H1_1*dH2_2 - c1*dH1_1*H2_2;
  T(1,2) = H2_1*dH2_2 - c1*dH2_1*H2_2;
  T(1,3) = c2*H1_1*H2_2;
  T(1,4) = c2*H2_1*H2_2;
  
  T(2,1) = c1*dH1_1*H1_2 - H1_1*dH1_2;
  T(2,2) = c1*dH2_1*H1_2 - H2_1*dH1_2;
  T(2,3) = -c2*H1_1*H1_2;
  T(2,4) = -c2*H2_1*H1_2;
  
  T(3,1) = d2*H1_1*H2_2;
  T(3,2) = d2*H2_1*H2_2;
  T(3,3) = H1_1*dH2_2 - d1*dH1_1*H2_2;
  T(3,4) = H2_1*dH2_2 - d1*dH2_1*H2_2;
  
  T(4,1) = -d2*H1_1*H1_2;
  T(4,2) = -d2*H2_1*H1_2;
  T(4,3) = d1*dH1_1*H1_2 - H1_1*dH1_2;
  T(4,4) = d1*dH2_1*H1_2 - H2_1*dH1_2;

  T *= pi * rk2 / (-4.*I);

  if (scaling)
  {
    cMatrix S1(4,4,fortranArray);
    cMatrix invS2(4,4,fortranArray);
    S1 = 0.0;
    invS2 = 0.0;
    
    S1(1,1) = S1(3,3) = exp(+2.0*I*rk1);
    S1(2,2) = S1(4,4) = 1.0;

    invS2(1,1) = invS2(3,3) = exp(-2.0*I*rk2);
    invS2(2,2) = invS2(4,4) = 1.0;

    T = multiply(invS2, T, S1);
    T = 10.0 * T; // For stability reasons, temporary fix.
  }

  return T;
  
}


/////////////////////////////////////////////////////////////////////////////
//
// Calculates field matrix.
// 
// To get the actual field one multiplies the field matrix with the 
// amplitude vector. The 6-component vector of field components has the 
// ordering [Ez Ephi Er Hz Hphi Hr]
//
/////////////////////////////////////////////////////////////////////////////

cMatrix field_matrix(const Complex& r, 
		     const Material& mat,
		     const Complex& kz_,
		     const int& circ_order,
		     const bool& scaling)
{

  const Complex kz = - kz_; // exp(-ikz) in CAMFR, versus exp(ikz) in Matlab
  
  // Set constants

  const Real    k0  = - 2*pi/global.lambda; // exp(-iwt) versus exp(iwt)
  const Complex ord = circ_order;

  const Complex n = mat.n();

  Complex k  = sqrt(k0*k0*n*n - kz*kz);
  // always put kt in the semicircle centered on the first quadrant
  if ((real(k)+imag(k))<0)
     {
    k = -k;
    // give a warning if k is close to this branch-cut
    if (abs(real(k)+imag(k))<0.01*abs(k))
      cout << "Warning: branch-cut in circMringutil.cpp: " << -k 
           << "became " << k << endl;
  }
  const Complex rk = k*r;

  const Complex c1 = -kz/k; 
  const Complex c4 = -k0/k; 
  const Complex c2 =  I*c4;
  const Complex c3 = -I*c1;
  const Complex d2 = -c2*n*n;
  const Complex d4 = -c4*n*n;

  // Calculate Bessel functions
  
  Complex H1, dH1dr;
  Complex H2, dH2dr;  
  Complex H1r;
  Complex H2r;

  if (rk != 0.0)
  {
    dH1dr = dH1(circ_order, rk, &H1, NULL, scaling);
    dH2dr = dH2(circ_order, rk, &H2, NULL, scaling);

    H1r = H1*ord/rk;
    H2r = H2*ord/rk;
  }
  else
  {
    Real J_, dJdr = dJ(circ_order, 0.0, &J_, NULL, scaling);
    H1    = J_;
    dH1dr = dJdr;    
    H2    = J_;
    dH2dr = dJdr;
    
    Real Jplus  = J(circ_order+1, 0.0, scaling);
    Real Jminus;
    if (circ_order >= 1)
      Jminus = J(circ_order-1, 0.0, scaling); 
    else
      Jminus = 0;
    H1r = 0.5 * (Jplus + Jminus);    
    H2r = H1r;
  }
  
  // Calculate field matrix

  cMatrix F(6,4,fortranArray);

  F(1,1) = H1;       F(1,2) = H2;       F(1,3) = 0.0;      F(1,4) = 0.0;
  F(2,1) = c1* H1r;  F(2,2) = c1* H2r;  F(2,3) = c2*dH1dr; F(2,4) = c2*dH2dr;
  F(3,1) = c3*dH1dr; F(3,2) = c3*dH2dr; F(3,3) = c4* H1r;  F(3,4) = c4* H2r;
  F(4,1) = 0.0;      F(4,2) = 0.0;      F(4,3) = H1;       F(4,4) = H2;
  F(5,1) = d2*dH1dr; F(5,2) = d2*dH2dr; F(5,3) = c1* H1r;  F(5,4) = c1* H2r;
  F(6,1) = d4* H1r;  F(6,2) = d4* H2r;  F(6,3) = c3*dH1dr; F(6,4) = c3*dH2dr;

  if (scaling == true)
  {
    cMatrix S(4,4,fortranArray);
    S = 0;

    S(1,1) = S(3,3) = exp(+2.0*I*rk);
    S(2,2) = S(4,4) = 1.0;

    F = multiply(F, S);
  }

  cMatrix Units(6,6,fortranArray);
  Units = 0;
  for (int i=1; i <= 3; i++) 
    Units(i, i) = 1.0; // electric field components
  
  const Complex mu0 = 4*pi*1e-07; 

  for (int i=4; i <= 6; i++)
    Units(i, i) = 1.0/(c*mu0); // magnetic field components

  F = multiply(Units, F);

  return F;
}
