
/////////////////////////////////////////////////////////////////////////////
//
// File:     circdisp.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990824
// Version:  1.5
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CIRCDISP_H
#define CIRCDISP_H

#include "../../defs.h"
#include "../../material.h"
#include "../../mode.h"
#include "../../math/calculus/calculus.h"
#include "../../math/bessel/bessel.h"

/////////////////////////////////////////////////////////////////////////////
//
// typedefs
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {kind_1, kind_2} Hankel;
typedef enum {guided, rad}    Guided_rad;



/////////////////////////////////////////////////////////////////////////////
//
// Function objects for dispersion relations for circular geometries
// with a homogeneously filled metal cylinder.
// Uses kr*rho as independent variable.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// TE modes
//
/////////////////////////////////////////////////////////////////////////////

class Circ_1_TE : public RealFunction
{
  public:

    Circ_1_TE(int n) : order(n) {}
    
    Real operator()(const Real& kr_rho)
      {counter++; return dJ(order, kr_rho);}
    
  protected:

    int order;
};



/////////////////////////////////////////////////////////////////////////////
//
// TM modes
//
/////////////////////////////////////////////////////////////////////////////

class Circ_1_TM : public RealFunction
{
  public:

    Circ_1_TM(int n) : order(n) {}

    Real operator()(const Real& kr_rho)
      {counter++; return J(order,kr_rho);}
    
  protected:

    int order;
};



/////////////////////////////////////////////////////////////////////////////
//
// Function objects for dispersion relations for circular geometries,
// using two materials and with support for complex coordinate stretching.
// Either open or closed systems and well-behaved for finding either
// guided modes or radiation/leaky modes.
// Uses kr2 as independent variable.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// Open system.
//
/////////////////////////////////////////////////////////////////////////////

class Circ_2_open : public ComplexFunction
{
  public:

    Circ_2_open(const Complex&  _r,
                const Material& _core,   const Material& _cladding,
                Real            _lambda, int             _order,
                Guided_rad      _type=guided)

      : r(_r), core(_core), cladding(_cladding),
        lambda(_lambda), order(_order), type(_type) {}

    Complex operator()(const Complex& kr2);

    vector<Complex> get_params() const;
    void            set_params(const vector<Complex>&);
    
  protected:

          Complex    r;         // core radius
          Material   core;
          Material   cladding;
          Real       lambda;    // wavelength
    const int        order;     // order of Bessel functions
    const Guided_rad type;      // variations
};



/////////////////////////////////////////////////////////////////////////////
//
// Closed system, always zero at kr1=0 and kr2=0.
//
/////////////////////////////////////////////////////////////////////////////

Complex scale(Complex* J, Complex* dJ, const Complex& z, Hankel hankel);

class Circ_2_closed : public ComplexFunction
{
  public:
  
    Circ_2_closed(const Complex&  _r,      const Complex&  _R,
                  const Material& _core,   const Material& _cladding,
                  Real            _lambda, int             _order,
                  Guided_rad      _type=guided,
                  Hankel          _hankel=kind_1,
                  Polarisation    _pol_0=TE,
                  bool            scale_always=false);

    void ab_TE_TM(Complex* a_TE, Complex* b_TE,
                  Complex* a_TM, Complex* b_TM,
                  const Complex& kr2) const;

    Complex operator()(const Complex& kr2);

    vector<Complex> get_params() const;
    void            set_params(const vector<Complex>&);
        
  protected:

          Complex      r;         // core radius
          Complex      R;         // metal cylinder radius
          Material     core;
          Material     cladding;
          Real         lambda;    // wavelength
    const int          order;     // order of Bessel functions
    const Guided_rad   type;      // variations
    const Hankel       hankel;    // kind of Hankel functions
          Polarisation pol_0;     // simplified for TE or TM order 0
    
    // exponential scaling of bessel functions in intermediate results
    
    bool scaling_co;    // core
    bool scaling_cl;    // cladding
    bool scaling_split; // ony scale in 'interesting' region  
};



/////////////////////////////////////////////////////////////////////////////
//
// Closed system, only zero at kr2=0 for true cut-off mode.
//
/////////////////////////////////////////////////////////////////////////////

class Circ_2_closed_cutoff : public Circ_2_closed
{
  public:
  
    Circ_2_closed_cutoff(const Complex&  r,      const Complex&  R,
                         const Material& core,   const Material& cladding,
                         Real            lambda, int             order,
                         Guided_rad      type=guided,
                         Hankel          hankel=kind_1,
                         Polarisation    pol_0=TE,
                         bool            scale_always=false)

      : Circ_2_closed(r,R,core,cladding,lambda,order,
                      rad,hankel,pol_0,scale_always) {}

    Complex operator()(const Complex& kr2);
};



/////////////////////////////////////////////////////////////////////////////
//
// Closed system, well behaved for finding radiation modes.
// Rotates results to be real-valued for guided modes in the case of
// exponential scaling in lossless structures.
//
/////////////////////////////////////////////////////////////////////////////

class Circ_2_closed_rad_lossless : public Circ_2_closed
{
  public:
  
    Circ_2_closed_rad_lossless
      (const Complex&  r,      const Complex&  R,
       const Material& core,   const Material& cladding,
       Real            lambda, int             order,
       Hankel          hankel=kind_1,
       Polarisation    pol_0=TE,
       bool            scale_always=false)

      : Circ_2_closed(r,R,core,cladding,lambda,order,
                      rad,hankel,pol_0,scale_always) {}

    Complex operator()(const Complex& kr2);
};


#endif
