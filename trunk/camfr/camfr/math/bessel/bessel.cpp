
/////////////////////////////////////////////////////////////////////////////
//
// File:     bessel.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990430
// Version:  1.02
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include "bessel.h"

using namespace std;

const Real eps = 1e-15;

/////////////////////////////////////////////////////////////////////////////
//
// BesselJ0, real argument x
//
//   low-level, not callable from the outside 
//
/////////////////////////////////////////////////////////////////////////////

extern "C" double F77NAME(dbesj0)(double&);

Real J0(Real x)
{
  return F77NAME(dbesj0)(x);
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselJ1, real argument x
//
//   low-level, not callable from the outside 
//
/////////////////////////////////////////////////////////////////////////////

extern "C" double F77NAME(dbesj1)(double&);

Real J1(Real x)
{
  return F77NAME(dbesj1)(x);
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselY0, real argument x
//
//   low-level, not callable from the outside 
//
/////////////////////////////////////////////////////////////////////////////

extern "C" double F77NAME(dbesy0)(double&);

Real Y0(Real x)
{
  return F77NAME(dbesy0)(x);
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselY1, real argument x
//
//   low-level, not callable from the outside 
//
/////////////////////////////////////////////////////////////////////////////

extern "C" double F77NAME(dbesy1)(double&);

Real Y1(Real x)
{
  return F77NAME(dbesy1)(x);
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselJ, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(dbesj)(double&, double&, int&, double*, int&);

Real J(Real n, Real x, bool scaled)
{
  if (n==0)
     return J0(x);
  
  if (n==1)
     return J1(x);
  
  int  orders=1;       // orders=1 : only compute this order
  Real result[1];
  int  underflows=0;
  
  F77NAME(dbesj)(x,n,orders,result,underflows);

  if (underflows)
  {
    std::ostringstream s;
    s << "Warning: " << underflows << " underflow(s) in J("
      << "("<< n << "," << x << ")";
    py_print(s.str());
  }

  return result[0];  
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselY, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(dbesy)(double&, double&, int&, double*);

Real Y(Real n, Real x, bool scaled)
{
  if (n==0)
     return Y0(x);
  
  if (n==1)
     return Y1(x);

  int  orders=1;       // orders=1 : only compute this order
  Real result[1];
  
  F77NAME(dbesy)(x,n,orders,result);

  return result[0]; 
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselH1, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

const Complex H1(Real n, Real x, bool scaled)
{  
  const Complex I(0,1);

  Complex factor(1,0);
  if (scaled) // make consistent with scaling in complex case
    factor = exp(-I*x);

  return (J(n,x,scaled) + I*Y(n,x,scaled)) * factor;  
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselH2, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

const Complex H2(Real n, Real x, bool scaled)
{
  const Complex I(0,1);
  
  Complex factor(1,0);
  if (scaled) // make consistent with scaling in complex case
    factor = exp(+I*x);
  
  return (J(n,x,scaled) - I*Y(n,x,scaled)) * factor;
}



/////////////////////////////////////////////////////////////////////////////
//
// errorcheck routine for complex Bessel functions
//
/////////////////////////////////////////////////////////////////////////////

void checkerror(int underflows, int ierr, const char* function,
                Real n, const Complex& z)
{
  if ( (ierr==0) && (underflows==0))
    return;

  std::ostringstream s;

  if (ierr==1)
  {
    s << "Error: bad input in " << function << "(" << n << "," << z << ")";
    py_error(s.str());
    exit(-1);
  }

  if ( (ierr==2) || (ierr==4) || (ierr==5) )
  {
    s << "Error " << ierr << ": overflow in " << function << "("
      << n << "," << z << "), try exp. scaling.";
    py_error(s.str());
  }
    
  if (ierr==3)
    s << "Warning: large arguments in " << function << "("
      << n << "," << z << ")";
  
  if (underflows)
    s << "Warning: " << underflows << " underflow(s) in " << function
      << "("<< n << "," << z << ")";

  py_print(s.str());
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselJ, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zbesj)(double&, double&, double&, int&, int&,
                               double*, double*, int&, int&);

const Complex J(Real n, const Complex& z, bool scaled)
{
  Real zr=z.real();
  Real zi=z.imag();

  if (abs(zi) < eps) // real routines are faster
     return J(n,zr,scaled);

  int  code=(scaled ? 2 : 1); // code=2   : exponential scaling
  int  orders=1;              // orders=1 : only compute this order
  Real yr[1], yi[1];
  int  underflows=0;
  int  ierr=0;

  F77NAME(zbesj)(zr,zi,n,code,orders,yr,yi,underflows,ierr);

  checkerror(underflows,ierr,"J",n,z);
  
  return Complex(yr[0],yi[0]); 
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselY, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zbesy)(double&, double&, double&, int&, int&,
                               double*, double*, int&, double*, double*, int&);

const Complex Y(Real n, const Complex& z, bool scaled)
{  
  Real zr=z.real();
  Real zi=z.imag(); 

  if ( (abs(zi) < eps) && (zr > 0) )// real routines are faster
     return Y(n,zr,scaled);

  int  code=(scaled ? 2 : 1); // code=2   : exponential scaling
  int  orders=1;              // orders=1 : only compute this order
  Real workr[1], worki[1];
  Real yr[1],    yi[1];
  int  underflows=0;
  int  ierr=0;

  F77NAME(zbesy)(zr,zi,n,code,orders,yr,yi,underflows,workr,worki,ierr);

  checkerror(underflows,ierr,"Y",n,z);
    
  return Complex(yr[0],yi[0]); 
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselH1, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zbesh)(double&, double&, double&, int&, int&, int&,
                               double*, double*, int&, int&);

const Complex H1(Real n, const Complex& z, bool scaled)
{
  Real zr=z.real();
  Real zi=z.imag();

  if ( (abs(zi) < eps) && (zr > 0) ) // real routines are faster
     return H1(n,zr,scaled);

  int  code=(scaled ? 2 : 1); // code=2   : exponential scaling
  int  type=1;                // H1
  int  orders=1;              // orders=1 : only compute this order
  Real yr[1], yi[1];
  int  underflows=0;
  int  ierr=0;

  F77NAME(zbesh)(zr,zi,n,code,type,orders,yr,yi,underflows,ierr);

  checkerror(underflows,ierr,"H1",n,z);
  
  return Complex(yr[0],yi[0]); 
}



/////////////////////////////////////////////////////////////////////////////
//
// BesselH2, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

const Complex H2(Real n, const Complex& z, bool scaled)
{
  Real zr=z.real();
  Real zi=z.imag();

  if ( (abs(zi) < eps) && (zr > 0) ) // real routines are faster
     return H2(n,zr,scaled);

  int  code=(scaled ? 2 : 1); // code=2   : exponential scaling
  int  type=2;                // H2
  int  orders=1;              // orders=1 : only compute this order
  Real yr[1], yi[1];
  int  underflows=0;
  int  ierr=0;

  F77NAME(zbesh)(zr,zi,n,code,type,orders,yr,yi,underflows,ierr);

  checkerror(underflows,ierr,"H2",n,z);
    
  return Complex(yr[0],yi[0]); 
}



/////////////////////////////////////////////////////////////////////////////
//
// Derivatives are calculated using the backward recursion formula:
//
//   x dZ(n,x) = - n Z(n,x) + x Z(n-1,x)
//
// The forward recursion formula is usually slower:
//
//   x dZ(n,x) =   n Z(n,x) - x Z(n+1,x)
//
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselJ, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

Real dJ(Real n, Real x, Real* Jn, Real* Jn_1, bool scaled)
{
  if (n==0)
  {
    Real J1=J(1,x,scaled);
    
    if (Jn)
      *Jn=J(0,x,scaled);

    if (Jn_1)
      *Jn_1=-J1; // J(0-1,x) = - J(1,x)
    
    return -J1; // dJ(0,x) = -J(1,x)
  }
  
  if (n==1) // specific routines are faster
  {
    Real J0=J(0,x,scaled);
    Real J1=J(1,x,scaled);

    if (Jn)
      *Jn=J1;

    if (Jn_1)
      *Jn_1=J0;

    if (abs(x) < eps)
      return 0.5;
    else
      return -J1/x + J0;
  }
  
  double startorder=n-1; // calculate orders n-1 and
  int    orders=2;       // .. and n
  int    underflows=0;
  Real   result[2];

  F77NAME(dbesj)(x,startorder,orders,result,underflows);

  if (underflows)
  {
    std::ostringstream s;
    s << "Warning: " << underflows << " underflow(s) in J("
      << "("<< n << "," << x << ")";
    py_error(s.str());
  }
  
  Real jn_1=result[0];
  Real   jn=result[1];
  
  if (Jn)
    *Jn=jn;

  if (Jn_1)
    *Jn_1=jn_1;

  if (abs(x) < eps)
    return 0;
  else
    return -jn*n/x + jn_1;  
}



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselY, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

Real dY(Real n, Real x, Real* Yn, Real* Yn_1, bool scaled)
{
  if (n==0)
  {
    Real Y1=Y(1,x,scaled);
    
    if (Yn)
      *Yn=Y(0,x,scaled);

    if (Yn_1)
      *Yn_1=-Y1; // Y(0-1,x) = - Y(1,x)
    
    return -Y1; // dY(0,x) = -Y(1,x)
  }
  
  if (n==1) // specific routines are faster
  {
    Real Y0=Y(0,x,scaled);
    Real Y1=Y(1,x,scaled);

    if (Yn)
      *Yn=Y1;

    if (Yn_1)
      *Yn_1=Y0;

    return -Y1/x + Y0;
  }
     
  double startorder=n-1; // calculate orders n-1 and
  int    orders=2;       // .. and n
  Real   result[2];

  F77NAME(dbesy)(x,startorder,orders,result);

  Real yn_1=result[0];
  Real   yn=result[1];
  
  if (Yn)
    *Yn=yn;

  if (Yn_1)
    *Yn_1=yn_1;
  
  return -yn*n/x + yn_1;   
}



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselH1, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////

const Complex dH1(Real n, Real x, Complex* Hn, Complex* Hn_1, bool scaled)
{
  Real dj, Jn, Jn_1;
  dj = dJ(n,x,&Jn,&Jn_1,scaled);

  Real dy, Yn, Yn_1;
  dy = dY(n,x,&Yn,&Yn_1,scaled);
  
  const Complex I(0,1);
  
  Complex factor(1,0);  
  if (scaled) // make consistent with scaling in complex case
    factor = exp(-I*x);
  
  if (Hn)
    *Hn = (Jn + I*Yn) * factor;

  if (Hn_1)
    *Hn_1 = (Jn_1 + I*Yn_1) * factor;

  return (dj + I*dy) * factor;  
}


  
/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselH2, real order n, real argument x
//
/////////////////////////////////////////////////////////////////////////////
  
const Complex dH2(Real n, Real x, Complex* Hn, Complex* Hn_1, bool scaled)
{
  Real dj, Jn, Jn_1;
  dj = dJ(n,x,&Jn,&Jn_1,scaled);

  Real dy, Yn, Yn_1;
  dy = dY(n,x,&Yn,&Yn_1,scaled);
  
  const Complex I(0,1);
  
  Complex factor(1,0);
  if (scaled) // Make consistent with scaling in complex case.
    factor = exp(+I*x);

  if (Hn)
    *Hn = (Jn - I*Yn) * factor;

  if (Hn_1)
    *Hn_1 = (Jn_1 - I*Yn_1) * factor;

  return (dj - I*dy) * factor; 
}



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselJ, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

const Complex dJ(Real n, const Complex& z, Complex* Jn, Complex* Jn_1,
                 bool scaled)
{
  if (n==0)
  {
    Complex J1=J(1,z,scaled);
    
    if (Jn)
      *Jn=J(0,z,scaled);

    if (Jn_1)
      *Jn_1=-J1; // J(0-1,z) = - J(1,z)
    
    return -J1; // dJ(0,z) = -J(1,z)
  }

  Real zr=z.real();
  Real zi=z.imag();

  if (abs(zi) < eps) // real routines are faster
  {
    Real dj, jn, jn_1;

    dj = dJ(n,zr,&jn,&jn_1,scaled);

    if (Jn)
      *Jn=jn;

    if (Jn_1)
      *Jn_1=jn_1;
    
    return dj;
  }
  
  int    code=(scaled ? 2 : 1); // code=2 : exponential scaling
  double startorder=n-1;        // calculate orders n-1 and
  int    orders=2;              // .. and n
  Real   yr[2], yi[2];
  int    underflows=0;
  int    ierr=0;

  F77NAME(zbesj)(zr,zi,startorder,code,orders,yr,yi,underflows,ierr);

  checkerror(underflows,ierr,"J",n,z);
  
  Complex jn_1(yr[0],yi[0]);
  Complex   jn(yr[1],yi[1]);

  if (Jn)
    *Jn=jn;

  if (Jn_1)
    *Jn_1=jn_1;

  if (abs(z) < eps)
    if (n==1)
      return 0.5;
    else
      return 0;
  else
    return -jn*n/z + jn_1;  
}



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselY, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

const Complex dY(Real n, const Complex& z, Complex* Yn, Complex* Yn_1,
                 bool scaled)
{
  if (n==0)
  {
    Complex Y1=Y(1,z,scaled);
    
    if (Yn)
      *Yn=Y(0,z,scaled);

    if (Yn_1)
      *Yn_1=-Y1; // Y(0-1,z) = - Y(1,z)
    
    return -Y1; // dY(0,z) = -Y(1,z)
  }

  Real zr=z.real();
  Real zi=z.imag();

  if ( (abs(zi) < eps) && (zr > 0) ) // real routines are faster
  {
    Real dy, yn, yn_1;

    dy = dY(n,zr,&yn,&yn_1,scaled);

    if (Yn)
      *Yn=yn;

    if (Yn_1)
      *Yn_1=yn_1;
    
    return dy;
  }
 
  int    code=(scaled ? 2 : 1); // code=2 : exponential scaling
  double startorder=n-1;        // calculate orders n-1 and
  int    orders=2;              // .. and n
  Real   workr[2], worki[2];
  Real   yr[2], yi[2];
  int    underflows=0;
  int    ierr=0;

  F77NAME(zbesy)(zr,zi,startorder,code,orders,yr,yi,
                 underflows,workr,worki,ierr);

  checkerror(underflows,ierr,"Y",n,z);
    
  Complex yn_1(yr[0],yi[0]);
  Complex   yn(yr[1],yi[1]);
  
  if (Yn)
    *Yn=yn;

  if (Yn_1)
    *Yn_1=yn_1;

  return -yn*n/z + yn_1;  
}


  
/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselH1, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

const Complex dH1(Real n, const Complex& z, Complex* Hn, Complex* Hn_1,
                  bool scaled)
{
  if (n==0)
  {
    Complex H1_1=H1(1,z,scaled);
    
    if (Hn)
      *Hn=H1(0,z,scaled);

    if (Hn_1)
      *Hn_1=-H1_1; // H1(0-1,z) = - H1(1,z)
    
    return -H1_1; // dH1(0,z) = -H1(1,z)
  }

  Real zr=z.real();
  Real zi=z.imag();

  if ( (abs(zi) < eps) && (zr > 0) ) // real routines are faster
     return dH1(n,zr,Hn,Hn_1,scaled);

  int    code=(scaled ? 2 : 1); // code=2 : exponential scaling
  int    type=1;                // H1
  double startorder=n-1;        // calculate orders n-1 and
  int    orders=2;              // .. and n
  Real   yr[2], yi[2];
  int    underflows=0;
  int    ierr=0;

  F77NAME(zbesh)(zr,zi,startorder,code,type,orders,yr,yi,underflows,ierr);
  
  checkerror(underflows,ierr,"H1",n,z);
  
  Complex H1n_1(yr[0],yi[0]);
  Complex   H1n(yr[1],yi[1]);

  if (Hn)
    *Hn=H1n;

  if (Hn_1)
    *Hn_1=H1n_1;

  return -H1n*n/z + H1n_1;  
}



/////////////////////////////////////////////////////////////////////////////
//
// derivative of BesselH2, real order n, complex argument z
//
/////////////////////////////////////////////////////////////////////////////

const Complex dH2(Real n, const Complex& z, Complex* Hn, Complex* Hn_1,
                  bool scaled)
{  
  if (n==0)
  {
    Complex H2_1=H2(1,z,scaled);
    
    if (Hn)
      *Hn=H2(0,z,scaled);

    if (Hn_1)
      *Hn_1=-H2_1; // H2(0-1,z) = - H2(1,z)
    
    return -H2_1; // dH2(0,z) = -H2(1,z)
  }

  Real zr=z.real();
  Real zi=z.imag();

  if ( (abs(zi) < eps) && (zr > 0) ) // real routines are faster
     return dH2(n,zr,Hn,Hn_1,scaled);

  int    code=(scaled ? 2 : 1); // code=2 : exponential scaling
  int    type=2;                // H2
  double startorder=n-1;        // calculate orders n-1 and
  int    orders=2;              // .. and n
  Real   yr[2], yi[2];
  int    underflows=0;
  int    ierr=0;

  F77NAME(zbesh)(zr,zi,startorder,code,type,orders,yr,yi,underflows,ierr);
  
  checkerror(underflows,ierr,"H2",n,z);
  
  Complex H2n_1(yr[0],yi[0]);
  Complex   H2n(yr[1],yi[1]);

  if (Hn)
    *Hn=H2n;

  if (Hn_1)
    *Hn_1=H2n_1;

  return -H2n*n/z + H2n_1;  
}


