
/////////////////////////////////////////////////////////////////////////////
//
// File:     circdisp.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990824
// Version:  1.5
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "circdisp.h"
#include "circ.h"

/////////////////////////////////////////////////////////////////////////////
//
// Circular geometry, 2 materials, open system.
// Support for complex coordinate stretching.
// Uses kr2 as independent variable.
// In this formulation with Hankel functions of second kind, guided modes
// correspond to neg imag(kr2), leaky modes to pos imag(kr2).
//
/////////////////////////////////////////////////////////////////////////////

Complex Circ_2_open::operator()(const Complex& kr2)
{
  counter++;

  // Set constants.

  const Real    eps   = 1e-15;
  
  const Complex epsr1 =     core.epsr();
  const Complex epsr2 = cladding.epsr();
  const Complex mur1  =     core.mur();
  const Complex mur2  = cladding.mur();
  
  const Real    k0    = 2*pi/lambda;
  const Complex beta2 = k0*k0*epsr2*mur2 - kr2*kr2;
  const Complex kr1   = signedsqrt(k0*k0*epsr1*mur1 - beta2, core);

  const bool  scaling = true;

  // Calculate bessel functions.
  
  Complex J1r, dJ1r = dJ(order, kr1*r, &J1r, NULL, scaling);
  Complex H2r, dH2r;
  
  if (abs(kr2*r) > eps)
    dH2r = dH2(order, kr2*r, &H2r, NULL, scaling);
  else
  {
    dH2r = 0; // limit kr2^2 dHr2 / H2r = 0
     H2r = 1; // if kr2 = 0
  }

  // Calculate dispersion relation.
  
  Complex T1, T2;
  
  if (type == rad) // Well behaved for radiation modes (no poles).
  {
    T1 = pow(r * kr1 * kr2, 2)

      * (   kr1 * epsr2 *  J1r * dH2r
          - kr2 * epsr1 * dJ1r *  H2r)

      * (   kr1 *  mur2 *  J1r * dH2r
          - kr2 *  mur1 * dJ1r *  H2r);

    T2 = -pow(Real(order) * J1r * H2r
              * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2;

    return T1 + T2;
  }
  
  if (type == guided) // Well behaved for guided modes (no exp. growth).
  {
    T1 = pow(r * kr1 * kr2, 2)

      * (   kr1 * epsr2 *  J1r * dH2r / H2r
          - kr2 * epsr1 * dJ1r)

      * (   kr1 *  mur2 *  J1r * dH2r / H2r
          - kr2 *  mur1 * dJ1r);

    T2 = -pow(Real(order) * J1r
              * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2;

    return T1 + T2;
  }

  cerr << "Invalid type for dispersion relation." << endl;
  exit (-1); 
}



/////////////////////////////////////////////////////////////////////////////
//
// get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Circ_2_open::get_params() const
{ 
  vector<Complex> params;

  params.push_back(r);
  params.push_back(core.n());
  params.push_back(core.mur());
  params.push_back(cladding.n());
  params.push_back(cladding.mur());
  params.push_back(lambda);

  return params;  
}


/////////////////////////////////////////////////////////////////////////////
//
// set_params
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2_open::set_params(const vector<Complex>& params)
{
  r        = params[0];
  core     = Material(params[1], params[2]);
  cladding = Material(params[3], params[4]);
  lambda   = real(params[5]);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2_closed::Circ_2_closed
//
/////////////////////////////////////////////////////////////////////////////

Circ_2_closed::Circ_2_closed(const Complex&  _r,    const Complex&  _R,
                             const Material& _core, const Material& _cladding,
                             Real            _lambda, int           _order,
                             Guided_rad      _type,
                             Hankel          _hankel,
                             Polarisation    _pol_0,
                             bool            scale_always)
  : r(_r), R(_R), core(_core), cladding(_cladding), lambda(_lambda),
    order(_order), type(_type), hankel(_hankel), pol_0(_pol_0)
{
  if (r==R)
  {
    cerr << "Use Circ_1_closed for uniform structures." << endl;
    exit (-1);
  }

  if ( (order == 0) && (pol_0 != TE) && (pol_0 != TM) )
  {
    cout << "Invalid polarisation for this order, changing to TE." << endl;
    pol_0 = TE;
  }
  
  // Split scaling only performs cladding scaling in 'interesting' region
  // (pos imag kr2 for H1, neg for H2).
  // This prevents a vanishing function in the 'uninteresting' region,
  // which could be mistaken for zeros.

  scaling_co = true;
  
  if (     scale_always
       || (real(r)/lambda > 25)
       || (real(R)/lambda > 25) ) // heuristic
  {
    scaling_cl    = true;
    scaling_split = true;
  }
  else // scaling (slower!) not needed
  {
    scaling_cl    = false;
    scaling_split = false;
  } 
}



/////////////////////////////////////////////////////////////////////////////
//
// Helper function for scaling purposes.
//
// On input J and dJ should contain scaled values according to
//    J_scal' = J(z).exp(-abs(z_imag))
// On output, these values will be scaled according to
//    J_scal  = J(z).exp(-I.z) ( for neg imag(z) )
//              J(z).exp(+I.z) ( for pos imag(z) )
//
// The function returns s given by (d)J/(d)H = (d)J_scal/(d)H_scal.exp(s),
// for scaled Hankel functions of kind 'hankel'.
//
/////////////////////////////////////////////////////////////////////////////

Complex scale(Complex* J, Complex* dJ, const Complex& z, Hankel hankel)
{ 
  if (imag(z) < 0.0)
  {
    *J  *= exp(-I*real(z));
    *dJ *= exp(-I*real(z));
    
    return (hankel == kind_1) ? 0.0 : 2.0*I*z;
  }
  else
  {
    *J  *= exp(+I*real(z));
    *dJ *= exp(+I*real(z));

    return (hankel == kind_1) ? -2.0*I*z : 0.0;
  }  
}



/////////////////////////////////////////////////////////////////////////////
//
// Calculates some constants used in dispersion relation.
// These calculations are stable for pos imag arguments when using Hankel
// functions of the first kind and stable for neg imag arguments for
// second kind.
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2_closed::ab_TE_TM(Complex* a_TE, Complex* b_TE,
                             Complex* a_TM, Complex* b_TM,
                             const Complex& kr2) const
{
  // Determine whether we should scale.
  
  bool scaling_here = scaling_cl;

  if (scaling_cl && scaling_split)
    if  (    ( (hankel == kind_1) && (imag(kr2) < 0) )
          || ( (hankel == kind_2) && (imag(kr2) > 0) ) )
        scaling_here = false;

  // Calculate bessel functions.


  Complex J2r, dJ2r, H2r, dH2r;
  Complex J2R, dJ2R, H2R, dH2R;

  // Calculate J(z).exp(-abs(z_imag)).
  
  dJ2r = dJ(order, kr2*r, &J2r, NULL, scaling_here);
  dJ2R = dJ(order, kr2*R, &J2R, NULL, scaling_here);

  // Calculate scaling factors. Note that scaling leads to
  // functions that are not holomorphic across imag(z) = 0.
  
  Complex s_r = 0.0;
  Complex s_R = 0.0;
  
  if (scaling_here)
  {
    s_r = scale(&J2r, &dJ2r, kr2*r, hankel);
    s_R = scale(&J2R, &dJ2R, kr2*R, hankel);      
  }

  const Complex s = exp(s_r - s_R);
  
  // Calculate H1(z).exp(-I.z) or H2(z).exp(+I.z).

  if (hankel == kind_1)
  {
    dH2r = dH1(order, kr2*r, &H2r, NULL, scaling_here);
    dH2R = dH1(order, kr2*R, &H2R, NULL, scaling_here);
  }
  else
  {
    dH2r = dH2(order, kr2*r, &H2r, NULL, scaling_here);
    dH2R = dH2(order, kr2*R, &H2R, NULL, scaling_here);
  }

  // Calculate results.

  if (a_TE && b_TE && a_TM && b_TM)
  {
    *a_TE = dJ2R * dH2r - dH2R * dJ2r * s;
    *b_TE = dJ2R *  H2r - dH2R *  J2r * s;
    
    *a_TM =  J2R * dH2r -  H2R * dJ2r * s;
    *b_TM =  J2R *  H2r -  H2R *  J2r * s; 
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Circular geometry, 2 materials, closed system.
// Support for complex coordinate stretching.
// Uses kr2 as independent variable.
// Always zero at kr1=0 and kr2=0.
// Useable far from kr2=0.
//
/////////////////////////////////////////////////////////////////////////////

Complex Circ_2_closed::operator()(const Complex& kr2_)
{
  counter++;

  // This function is mathematically even in kr2. However, it is only
  // numerically stable in one half plane.

  Complex kr2 = kr2_;

  if (   ( (hankel == kind_1) && (imag(kr2) < 0) )
       ||( (hankel == kind_2) && (imag(kr2) > 0) ) )
    kr2 = -kr2;
  
  // Set constants.

  const Real    eps   = 1e-15;
  
  const Complex epsr1 =     core.epsr();
  const Complex epsr2 = cladding.epsr();
  const Complex mur1  =     core.mur();
  const Complex mur2  = cladding.mur();
  
  const Real    k0    = 2*pi/lambda;
  const Complex beta2 = k0*k0*epsr2*mur2 - kr2*kr2;
  const Complex kr1   = signedsqrt(k0*k0*epsr1*mur1 - beta2, core);

  // Limit for kr2=0.

  if (abs(kr2) < eps)
    return 0.0;
  
  // Calculate J(z).exp(-abs(z_imag))  (scaling factors cancel here).

  Complex J1r, dJ1r = dJ(order, kr1*r, &J1r, NULL, scaling_co);

  // Calculate dispersion relation.

  Complex a_TE, b_TE, a_TM, b_TM;
  ab_TE_TM(&a_TE, &b_TE, &a_TM, &b_TM, kr2);

  Complex T1, T2;
  
  if (order == 0)
  {
    if ( (pol_0 == TE) && (type == guided) )

      return I * (   kr1 *  mur2 *  J1r * a_TE / b_TE
                   - kr2 *  mur1 * dJ1r ) ;
    
    if ( (pol_0 == TE) && (type == rad) )

      return I * (   kr1 *  mur2 *  J1r * a_TE
                   - kr2 *  mur1 * dJ1r * b_TE );

    if ( (pol_0 == TM) && (type == guided) )

      return I * (   kr1 * epsr2 *  J1r * a_TM / b_TM
                   - kr2 * epsr1 * dJ1r );

    if ( (pol_0 == TM) && (type == rad) )

      return I * (   kr1 * epsr2 *  J1r * a_TM
                   - kr2 * epsr1 * dJ1r * b_TM );
  }
  else // order != 0
  {
    if (type == guided)
    {
      T1 = pow(r * kr1 * kr2, 2)

        * (   kr1 *  mur2 *  J1r * a_TE / b_TE
            - kr2 *  mur1 * dJ1r)

        * (   kr1 * epsr2 *  J1r * a_TM / b_TM
            - kr2 * epsr1 * dJ1r);

      T2 = -pow(Real(order) * J1r
                * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2;

      return T1 + T2;
    }
    
    if (type == rad)
    {
      T1 = pow(r * kr1 * kr2, 2)

        * (   kr1 *  mur2 *  J1r * a_TE
            - kr2 *  mur1 * dJ1r * b_TE )

        * (   kr1 * epsr2 *  J1r * a_TM
            - kr2 * epsr1 * dJ1r * b_TM );

      T2 = -pow(Real(order) * J1r
                * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2 * b_TE * b_TM;

      return T1 + T2;
    }
  }

  cerr << "Invalid type for dispersion relation." << endl;
  exit (-1);  
}



/////////////////////////////////////////////////////////////////////////////
//
// Circular geometry, 2 materials, closed system.
// Support for complex coordinate stretching.
// Uses kr2 as independent variable.
// Only zero at kr2=0 for true cut-off mode.
// Not useable far from kr2=0.
//
/////////////////////////////////////////////////////////////////////////////

Complex Circ_2_closed_cutoff::operator()(const Complex& kr2)
{
  counter++;

  // This function is mathematically even in kr2. However, it is only
  // numerically stable in one half plane.

  if (   ( (hankel == kind_1) && (imag(kr2) < 0) )
       ||( (hankel == kind_2) && (imag(kr2) > 0) ) )
    const_cast<Complex&>(kr2) = -kr2;
  
  // Set constants.

  const Real    eps   = 1e-15;
  
  const Complex epsr1 =     core.epsr();
  const Complex epsr2 = cladding.epsr();
  const Complex mur1  =     core.mur();
  const Complex mur2  = cladding.mur();
  
  const Real    k0    = 2*pi/lambda;
  const Complex beta2 = k0*k0*epsr2*mur2 - kr2*kr2;
  const Complex kr1   = signedsqrt(k0*k0*epsr1*mur1 - beta2, core);

  // Limit for kr2=0.

  if (abs(kr2) < eps)
    return 999; // inf
  
  // Calculate J(z).exp(-abs(z_imag))  (scaling factors cancel here).

  Complex J1r, dJ1r = dJ(order, kr1*r, &J1r, NULL, scaling_co);

  // Calculate dispersion relation.

  Complex a_TE, b_TE, a_TM, b_TM;
  ab_TE_TM(&a_TE, &b_TE, &a_TM, &b_TM, kr2);

  Complex T1, T2;
  
  if (order == 0)
  {
    if ( (pol_0 == TE) && (type == guided) )

      return I * (   1./kr2 *  mur2 *  J1r * a_TE / b_TE
                   - 1./kr1 *  mur1 * dJ1r ) ;
    
    if ( (pol_0 == TE) && (type == rad) )

      return I * (   1./kr2 *  mur2 *  J1r * a_TE
                   - 1./kr1 *  mur1 * dJ1r * b_TE );

    if ( (pol_0 == TM) && (type == guided) )

      return I * (   1./kr2 * epsr2 *  J1r * a_TM / b_TM
                   - 1./kr1 * epsr1 * dJ1r );

    if ( (pol_0 == TM) && (type == rad) )

      return I * (   1./kr2 * epsr2 *  J1r * a_TM
                   - 1./kr1 * epsr1 * dJ1r * b_TM );
  }
  else // order != 0
  {
    if (type == guided)
    {
      T1 = pow(r, 2)

        * (   1./kr2 *  mur2 *  J1r * a_TE / b_TE
            - 1./kr1 *  mur1 * dJ1r)

        * (   1./kr2 * epsr2 *  J1r * a_TM / b_TM
            - 1./kr1 * epsr1 * dJ1r);

      T2 = -pow(Real(order) * J1r
                * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2
                / pow(kr1 * kr2, 4);

      return T1 + T2;
    }
    
    if (type == rad)
    {
      T1 = pow(r, 2)

        * (   1./kr2 *  mur2 *  J1r * a_TE
            - 1./kr1 *  mur1 * dJ1r * b_TE )

        * (   1./kr2 * epsr2 *  J1r * a_TM
            - 1./kr1 * epsr1 * dJ1r * b_TM );

      T2 = -pow(Real(order) * J1r
                * k0 * (epsr1*mur1-epsr2*mur2), 2) * beta2
                * b_TE * b_TM / pow(kr1 * kr2, 4);

      return T1 + T2;
    }
  }

  cerr << "Invalid type for dispersion relation." << endl;
  exit (-1);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circular geometry, 2 materials, closed system.
// Support for complex coordinate stretching.
// Uses kr2 as independent variable.
// Well-behaved for finding radiation modes (no poles).
// Rotates results to be real-valued on axes for lossless structures.
//
/////////////////////////////////////////////////////////////////////////////

Complex Circ_2_closed_rad_lossless::operator()(const Complex& kr2)
{ 
  Complex result = Circ_2_closed::operator()(kr2);

  Complex sJ = 0.0;
  Complex sH = 0.0;
  
  if (scaling_cl)
  {
    sJ = (imag(kr2*R) < 0.0) ? I*kr2*R : -I*kr2*R;
    sH = (hankel == kind_1)  ? I*kr2*r : -I*kr2*r;
  }
  
  if (order == 0)
    return result*(exp(sJ+sH)); 
  else
    return result*(exp(2.0*(sJ+sH)));
}



/////////////////////////////////////////////////////////////////////////////
//
// get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Circ_2_closed::get_params() const
{   
  vector<Complex> params;

  params.push_back(r);
  params.push_back(R);
  params.push_back(core.n());
  params.push_back(core.mur());
  params.push_back(cladding.n());
  params.push_back(cladding.mur());
  params.push_back(lambda);

  return params;  
}



/////////////////////////////////////////////////////////////////////////////
//
// set_params
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2_closed::set_params(const vector<Complex>& params)
{
  r        = params[0];
  R        = params[1];
  core     = Material(params[2], params[3]);
  cladding = Material(params[4], params[5]);
  lambda   = real(params[6]);
}
