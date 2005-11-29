
/////////////////////////////////////////////////////////////////////////////
//
// File:     defs.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19980824
// Version:  1.0
//
// Copyright (C) 1998 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "defs.h"

Global global={0,0,TE,0,track,normal,100,1,0.01,100,100,Complex(1,1),false,
               20,1e-14,true,1e-12,identical,GEV,lapack,true,true,false,
               0.0,1.2,false,false,false,true};



/////////////////////////////////////////////////////////////////////////////
//
// Workaround for inferior compilers.
//
/////////////////////////////////////////////////////////////////////////////

#ifdef NO_CXX_ABS
Real abs(Real x) {return x >= 0.0 ? x : -x;}
#endif



/////////////////////////////////////////////////////////////////////////////
//
// operator<< for Polarisation
//
/////////////////////////////////////////////////////////////////////////////

const std::string Pol_string[] 
  = {"unknown", "TEM", "TE", "TM", "HE", "EH", "TE_TM"};

std::ostream& operator<< (std::ostream& s, const Polarisation& pol)
{
  return s << Pol_string[pol];
};



/////////////////////////////////////////////////////////////////////////////
//
// Python print functions.
//
/////////////////////////////////////////////////////////////////////////////

void py_print(const std::string& s) {PySys_WriteStdout("%s\n",s.c_str());}
void py_error(const std::string& s) {PySys_WriteStderr("%s\n",s.c_str());}



/////////////////////////////////////////////////////////////////////////////
//
// machine_eps
//
//   smallest number x for which 1+x != 1
//   routines from machar.c (see math/fortran/slatec)
//
/////////////////////////////////////////////////////////////////////////////

Real machine_eps()
{
  static Real calc_machine_eps = 0;
  
  if (calc_machine_eps) // already calculated
    return calc_machine_eps;

  //Real x=1, last_x;
  //
  //while (1+x != 1)
  //{
  //  last_x = x;
  //  x /= 2;
  //}
  //
  //calc_machine_eps = last_x;

  int i,iz,j,k;
  int mx,itmp,nxres;
  int ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd;
  Real a,b,beta,betain,one,y,z,zero;
  Real betah,t,tmp,tmpa,tmp1,two,epsneg;

  irnd = 1;
  one = Real(irnd);
  two = one + one;
  a = two;
  b = a;
  zero = 0.0e0;
  
  // Determine ibeta,beta ala malcolm.

  tmp = ((a+one)-a)-one;

  while (tmp == zero)
  {
    a = a+a;
    tmp = a+one;
    tmp1 = tmp-a;
    tmp = tmp1-one;
  }

  tmp = a+b;
  itmp = (int)(tmp-a);
  while (itmp == 0)
  {
    b = b+b;
    tmp = a+b;
    itmp = (int)(tmp-a);
  }

  ibeta = itmp;
  beta = Real(ibeta);
  
  // Determine irnd, it.

  it = 0;
  b = one;
  tmp = ((b+one)-b)-one;

  while (tmp == zero)
  {
    it = it+1;
    b = b*beta;
    tmp = b+one;
    tmp1 = tmp-b;
    tmp = tmp1-one;
  }

  irnd = 0;
  betah = beta/two;
  tmp = a+betah;
  tmp1 = tmp-a;
  if (tmp1 != zero) irnd = 1;
  tmpa = a+beta;
  tmp = tmpa+betah;
  if ((irnd == 0) && (tmp-tmpa != zero)) irnd = 2;
  
  // Determine negep, epsneg.

  negep = it + 3;
  betain = one / beta;
  a = one;

  for (i = 1; i<=negep; i++)
    a = a * betain;

  b = a;
  tmp = (one-a);
  tmp = tmp-one;

  while (tmp == zero)
  {
    a = a*beta;
    negep = negep-1;
    tmp1 = one-a;
    tmp = tmp1-one;
  }

  negep = -negep;
  epsneg = a;
  
  // Determine machep, eps.

  machep = -it - 3;
  a = b;
  tmp = one+a;

  while (tmp-one == zero)
  {
    a = a*beta;
    machep = machep+1;
    tmp = one+a;
  }

  calc_machine_eps = a;

  return calc_machine_eps;
}



/////////////////////////////////////////////////////////////////////////////
//
// out_of_memory
//   error handler
//
/////////////////////////////////////////////////////////////////////////////

void out_of_memory()
{
  py_error("Fatal error: Out of memory.");
  exit(-1);
}



/////////////////////////////////////////////////////////////////////////////
//
// pick_sign_k
//
/////////////////////////////////////////////////////////////////////////////

void pick_sign_k(Complex* k)
{
  // Lossy only.

  if (imag(*k) > 0)
    *k = - *k;

  if (abs(imag(*k)) < 1e-12)
    if (real(*k) < 0)
      *k = - *k;

  return;

  // Old style.

  if (real(*k) < 0)
    *k = - *k;

  if (abs(real(*k)) < 1e-8)
    if (imag(*k) > 0)
      *k = - *k;
}
