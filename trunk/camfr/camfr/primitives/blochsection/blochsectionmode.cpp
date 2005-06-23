
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsectionmode.cpp
// Author:   Peter.Bienstman@UGent.be
// Date:     20050518
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "../section/section.h"
#include "blochsectionmode.h"
#include "blochsectionoverlap.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

/////////////////////////////////////////////////////////////////////////////
//
// BlochSectionMode::BlochSectionMode
//
/////////////////////////////////////////////////////////////////////////////

BlochSectionMode::BlochSectionMode(Polarisation pol, const Complex& kz, 
                                   BlochSectionImpl* geom_,
                                   const cVector& Ex_, const cVector& Ey_,
                                   const cVector& Hx_, const cVector& Hy_)
: Mode(pol, kz, -kz), geom(geom_),
  Ex(fortranArray), Ey(fortranArray),
  Hx(fortranArray), Hy(fortranArray)
{
  Ex.resize(Ex_.shape()); Ex = Ex_;
  Ey.resize(Ey_.shape()); Ey = Ey_;
  Hx.resize(Hx_.shape()); Hx = Hx_;
  Hy.resize(Hy_.shape()); Hy = Hy_;
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSectionMode::field
//
/////////////////////////////////////////////////////////////////////////////

Field BlochSectionMode::field(const Coord& coord) const 
{
  Field f;

  Coord c = coord;
  c.c1 += I*global_section.left_PML;
  c.c2 += I*global_slab.lower_PML;

  const int M = global_blochsection.Mx;
  const int N = global_blochsection.My;

  const Complex W = get_geom()->get_width();
  const Complex H = get_geom()->get_height();

  const Complex alpha0 = global_blochsection.alpha0;
  const Complex  beta0 = global_blochsection.beta0;

  for (Real m=-M; m<=M; m+=1.0)
    for (Real n=-N; n<=N; n+=1.0)
    {
      int i1 = int((m+M+1) + (n+N)*(2*M+1));

      Complex alpha = alpha0 + m*2.*pi/W;
      Complex  beta =  beta0 + n*2.*pi/H;
      Complex expon = exp(I*(alpha*c.c1 + beta*c.c2));

      f.E1 += Ex(i1)*expon;
      f.E2 += Ey(i1)*expon;
      f.H1 += Hx(i1)*expon;
      f.H2 += Hy(i1)*expon;
    } 

    f.Ez = f.Hz = 0.0; // TMP.

    return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSectionMode::normalise
//
//  Note: because of complex conjugate, these modes cannot be normalised
//  e.g. when overlap(this,this) is complex. Also, complex modes have a zero
//  norm.
//
/////////////////////////////////////////////////////////////////////////////

void BlochSectionMode::normalise() 
{
  Complex norm = sqrt(overlap(this,this));
  
  if (abs(norm) < 1e-5)
  {
    //py_print("Warning: mode close to cutoff.");
    norm = 1.0;
  }
  
  Ex /= norm;
  Ey /= norm;
  Hx /= norm;
  Hy /= norm;

  return;

  std::cout << n_eff() << std::endl;
  std::cout << Ex << std::endl;
  std::cout << Ey << std::endl;
  std::cout << Hx << std::endl;
  std::cout << Hy << std::endl;  
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformBlochSectionMode::get_kx
//
/////////////////////////////////////////////////////////////////////////////

Complex UniformBlochSectionMode::get_kx() const
{
  const Complex W = get_geom()->get_width();
  const Complex alpha0 = global_blochsection.alpha0;
  return alpha0 + M*2.*pi/W;
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformBlochSectionMode::get_ky
//
/////////////////////////////////////////////////////////////////////////////

Complex UniformBlochSectionMode::get_ky() const
{
  const Complex  beta0 = global_blochsection.beta0;
  const Complex H = get_geom()->get_height();
  return beta0 + N*2.*pi/H;
}

