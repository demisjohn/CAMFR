
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
#include "circ_M_util.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

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
  
  const Complex k0    = 2*pi/lambda;
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
  lambda   = params[5];
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2_closed::Circ_2_closed
//
/////////////////////////////////////////////////////////////////////////////

Circ_2_closed::Circ_2_closed(const Complex&  _r,    const Complex&  _R,
                             const Material& _core, const Material& _cladding,
                             const Complex&  _lambda, int           _order,
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
       || (real(r)/real(lambda) > 25)
       || (real(R)/real(lambda) > 25) ) // heuristic
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
  
  const Complex k0    = 2*pi/lambda;
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
  
  const Complex k0    = 2*pi/lambda;
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
  lambda   = params[6];
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M_closed::Circ_M_closed
//
/////////////////////////////////////////////////////////////////////////////

Circ_M_closed::Circ_M_closed(unsigned int _M, 
			     const vector<Complex>&  _r,    
                             const vector<Material*>& _m, 
                             const Complex& _lambda, int _order,
                             Polarisation    _pol_0)
  : M(_M), radius(_r), material(_m),lambda(_lambda),
    order(_order), pol_0(_pol_0)
{
  if (M<3)
  {
    cout << "Error: use Circ_2 or Circ_1 for M<3" << endl;
    exit(-1);
  }

  if ( (order == 0) && (pol_0 != TE) && (pol_0 != TM) )
  {
    cout << "Invalid polarisation for this order, changing to TE." << endl;
    pol_0 = TE;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Circular geometry, M rings, closed system.
// Uses kz as independent variable.
// disp = 0 when boundary conditions are satisfied (Ez = Ephi = 0 at metal wall)
//
/////////////////////////////////////////////////////////////////////////////

Complex Circ_M_closed::operator()(const Complex& kt)
{
  counter++;

  const bool scaling = true; 

  const Complex nlast   = material[M-1]->n();
  const Complex murlast = material[M-1]->mur();
  const Complex k0      = 2*pi/global.lambda;
  const Complex klast   = k0*nlast*sqrt(murlast);

  Complex kz = sqrt(klast*klast - kt*kt);
  // always put kz in 4th or 1st quadrant
  if (real(kz)<0)
    kz = -kz;
  // also if kz is on imaginary axis make its imaginary part negative
  // (avoid noise problems)
  if (abs(real(kz))<1e-12)
    if (imag(kz)>0)
      kz = -kz;

  cMatrix Ttotal(4,4,fortranArray); 
  cMatrix F(6,4,fortranArray);
  cMatrix Mcore(4,2,fortranArray);
  cMatrix M_Etang(2,6,fortranArray);

  cMatrix Mresult(2,2,fortranArray);

  Ttotal = total_transfer_matrix(kz, scaling);

  F = field_matrix(radius[M-1], *material[M-1], kz, order, scaling);

  // Mcore says that at r = 0 the coefficients of the Y Bessel functions 
  // are zero

  Mcore(1,1) = 1.0;  Mcore(1,2) = 0.0;
  Mcore(2,1) = 1.0;  Mcore(2,2) = 0.0;
  Mcore(3,1) = 0.0;  Mcore(3,2) = 1.0;
  Mcore(4,1) = 0.0;  Mcore(4,2) = 1.0;

  M_Etang = 0.0;
  M_Etang(1,1) = 1.0; // select Ez   component
  M_Etang(2,2) = 1.0; // select Ephi component

  Mresult = multiply(M_Etang, F,  Ttotal, Mcore);

  if (order == 0) // temporary fix; for order == 0 
  // one should use 2x2 instead of 4x4 matrices
    if (pol_0 == TM)
      return Mresult(1,1);
    else
      return Mresult(2,2);
  
  return Mresult(1,1)*Mresult(2,2) - Mresult(1,2)*Mresult(2,1);
}


/////////////////////////////////////////////////////////////////////////////
//
// Calculates total transfer matrix for Circular geometry, M rings, 
// closed system.
// Uses kz as independent variable.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix Circ_M_closed::total_transfer_matrix(const Complex& kz, bool scaling)
{
  cMatrix T(4,4,fortranArray);

  T = 0.0;
  for (int i=1; i<=4; i++)
    T(i,i) = 1.0;

  for (unsigned int i=0; i<M-1; i++)
    T = multiply(transfer_matrix(radius[i], *material[i], *material[i+1], 
                                 kz, order, scaling), T);  
  return T;
}



/////////////////////////////////////////////////////////////////////////////
//
// get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Circ_M_closed::get_params() const
{   
  vector<Complex> params;

  params.push_back(lambda);
  for (unsigned int i=0; i<radius.size(); i++)
  {
    params.push_back(radius[i]);
    params.push_back(material[i]->n());
    params.push_back(material[i]->mur());
  }

  return params;  
}



/////////////////////////////////////////////////////////////////////////////
//
// set_params
//
/////////////////////////////////////////////////////////////////////////////

void Circ_M_closed::set_params(const vector<Complex>& params)
{
  lambda = params[0];

  unsigned int index = 0;

  for (unsigned int i=1; i<params.size(); i=i+3)
  {
    radius[index] = params[i];
    *material[index] = Material(params[i+1], params[i+2]);
    index++;
  }
}

