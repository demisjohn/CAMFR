
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
// BlochSection2D_Mode::BlochSection2D_Mode
//
/////////////////////////////////////////////////////////////////////////////

BlochSection2D_Mode::BlochSection2D_Mode
  (Polarisation pol, const Complex& kz, BlochSectionImpl* geom,
   const cVector& Ex_, const cVector& Ey_, 
   const cVector& Hx_, const cVector& Hy_)
    : BlochSectionMode(pol, kz, geom), 
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
// BlochSection2D_Mode::field
//
/////////////////////////////////////////////////////////////////////////////

Field BlochSection2D_Mode::field(const Coord& coord) const 
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

      Complex alpha  = alpha0 + m*2.*pi/W;
      Complex  beta  =  beta0 + n*2.*pi/H;
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
// BlochSection2D_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D_Mode::normalise() 
{
  Complex norm = sqrt(overlap(this, this));

  Ex /= norm;
  Ey /= norm;
  Hx /= norm;
  Hy /= norm;

  std::cout << kz << std::endl;
  std::cout << Ex << std::endl;;
  std::cout << Ey << std::endl;;
  std::cout << Hx << std::endl;;
  std::cout << Hy << std::endl;;
  std::cout << std::endl;
}


#if 0

/////////////////////////////////////////////////////////////////////////////
//
// BlochSection1D_Mode::BlochSection1D_Mode
//
/////////////////////////////////////////////////////////////////////////////

BlochSection1D_Mode::BlochSection1D_Mode
  (Polarisation pol, const Complex& kz, SlabMode* m_, BlochSection1D* geom) 
    : BlochSectionMode(pol, kz, geom), m(m_)
{
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection1D_Mode::field
//
/////////////////////////////////////////////////////////////////////////////

Field BlochSection1D_Mode::field(const Coord& coord) const 
{
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection1D_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection1D_Mode::normalise() 
{
}

#endif
