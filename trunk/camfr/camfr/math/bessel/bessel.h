
/////////////////////////////////////////////////////////////////////////////
//
// File:     bessel.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990129
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef BESSEL_H
#define BESSEL_H

#include <complex>
using namespace std;

typedef double          Real;
typedef complex<double> Complex;


/////////////////////////////////////////////////////////////////////////////
//
// Interface to the SLATEC Fortran routines computing Bessel functions
// (see www.netlib.org)
//
// Note on execution times: For complex arguments, Y takes twice as long
// as H and H takes takes twice as long as J.
// Complex arguments can take twice as long as real ones.
//
// When calculating derivatives dZn, there are two optional arguments that
// let you retrieve Zn and Zn-1 from the derivative calculation.
// This can avoid duplication of calculations.
//
// When scaled is set to true, the functions return
//
//      J(z).exp(-abs(z_imag))
//      Y(z).exp(-abs(z_imag))
//     H1(z).exp(-j.z)
//     H2(z).exp( j.z)
//
// and the same factors for the derivatives
//
/////////////////////////////////////////////////////////////////////////////

      Real     J(Real n, Real x, bool scaled=false);
      Real     Y(Real n, Real x, bool scaled=false);
const Complex H1(Real n, Real x, bool scaled=false);
const Complex H2(Real n, Real x, bool scaled=false);

const Complex  J(Real n, const Complex& z, bool scaled=false);
const Complex  Y(Real n, const Complex& z, bool scaled=false);
const Complex H1(Real n, const Complex& z, bool scaled=false);
const Complex H2(Real n, const Complex& z, bool scaled=false);

      Real     dJ(Real n, Real x, Real*    Jn=0, Real*    Jn_1=0,
                  bool scaled=false);
      Real     dY(Real n, Real x, Real*    Yn=0, Real*    Yn_1=0,
                  bool scaled=false);
const Complex dH1(Real n, Real x, Complex* Hn=0, Complex* Hn_1=0,
                  bool scaled=false);
const Complex dH2(Real n, Real x, Complex* Hn=0, Complex* Hn_1=0,
                  bool scaled=false);

const Complex  dJ(Real n, const Complex& z, Complex* Jn=0, Complex* Jn_1=0,
                  bool scaled=false);
const Complex  dY(Real n, const Complex& z, Complex* Yn=0, Complex* Yn_1=0,
                  bool scaled=false);
const Complex dH1(Real n, const Complex& z, Complex* Hn=0, Complex* Hn_1=0,
                  bool scaled=false);
const Complex dH2(Real n, const Complex& z, Complex* Hn=0, Complex* Hn_1=0,
                  bool scaled=false);



#endif




