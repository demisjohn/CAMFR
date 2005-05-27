
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsectionoverlap.cpp
// Author:   Peter.Bienstman@UGent.be
// Date:     20050518
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "blochsectionoverlap.h"
#include "blochsectionmode.h"
#include "../section/section.h"
#include "../../math/calculus/quadrature/patterson_quad.h"

using std::vector;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// Note: overlap integrals here are Int(E x H*) dz!
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// Overlap function for a given x. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_x_bloch : public ComplexFunction
{
  public:

    Overlap_x_bloch(const BlochSectionMode* m1_,const BlochSectionMode* m2_,
                    const Complex& x_)
      : m1(m1_), m2(m2_), x(x_) {}

    Complex operator()(const Complex& y)
    {
      counter++;

      Field f1 = m1->field(Coord(x,y,0));
      Field f2 = m2->field(Coord(x,y,0));

      return f1.E1 * conj(f2.H2) - f1.E2 * conj(f2.H1);
    }

  protected:

    const BlochSectionMode* m1;
    const BlochSectionMode* m2;
    Complex x;
};



/////////////////////////////////////////////////////////////////////////////
//
// Overlap function integrated out over y. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_bloch : public ComplexFunction
{
  public:

    Overlap_bloch(const BlochSectionMode* m1_, const BlochSectionMode* m2_) 
      : m1(m1_), m2(m2_) {}

    Complex operator()(const Complex& x)
    { 
      counter++;

      // Do numeric integration.

      Overlap_x_bloch f(m1, m2, x);

      Wrap_real_to_real f_r(f);
      Wrap_real_to_imag f_i(f);

      BlochSectionImpl* medium_I  
        = dynamic_cast<BlochSectionImpl*>(m1->get_geom());
      
      Complex H = medium_I->get_height();

      return     patterson_quad(f_r, 0.0, real(H), 1e-2, 4)
             + I*patterson_quad(f_i, 0.0, real(H), 1e-2, 4);
    }

  protected:

    const BlochSectionMode* m1;
    const BlochSectionMode* m2;
};



/////////////////////////////////////////////////////////////////////////////
//
// overlap_numeric_bloch
//
//   Verification function which calculates the overlap integrals 
//   numerically. Slower and less acurate.
//   TODO: check PML in z direction.
//   To be removed after the analytical overlap calculation matures.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap_numeric_bloch(const BlochSectionMode* mode_I,
                              const BlochSectionMode* mode_II)
{
  BlochSectionImpl* medium_I  
    = dynamic_cast<BlochSectionImpl*>(mode_I ->get_geom());

  Overlap_bloch f(mode_I, mode_II);

  Wrap_real_to_real f_r(f);
  Wrap_real_to_imag f_i(f);

  Complex W = medium_I->get_width();

  return     patterson_quad(f_r, 0.0, real(W), 1e-2, 4)
         + I*patterson_quad(f_i, 0.0, real(W), 1e-2, 4);
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const BlochSectionMode* sec_I_mode,
                const BlochSectionMode* sec_II_mode)
{ 
  const int M = global_blochsection.Mx;
  const int N = global_blochsection.My;

  const Complex alpha0 = global_blochsection.alpha0;
  const Complex  beta0 = global_blochsection.beta0;

  Complex W = sec_I_mode->get_geom()->get_width();
  Complex H = sec_I_mode->get_geom()->get_height();

  if (global_section.section_solver == L_anis)
  {
    W = real(W);
    H = real(H);
  }

  Complex result = 0.0;
  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);

      const Complex Ex = sec_I_mode->Ex(i1);
      const Complex Ey = sec_I_mode->Ey(i1);

      const Complex Hx = sec_II_mode->Hx(i1);
      const Complex Hy = sec_II_mode->Hy(i1);

      result += Ex*conj(Hy) - Ey*conj(Hx);
    }

  return result*W*H;

  Complex numeric = overlap_numeric_bloch(sec_I_mode, sec_II_mode);

  std::cout << "overlap " << numeric << " " << result*W*H
             << " " << numeric-result*W*H << std::endl;
  
  if (abs(numeric - result*W*H) > 1e-4)
    std::cout << "XXXX" << std::endl;
  
  return result*W*H;
}
