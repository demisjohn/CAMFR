
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectionoverlap.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020612
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "sectionoverlap.h"
#include "sectionmode.h"
#include "../../math/calculus/quadrature/patterson_quad.h"

/////////////////////////////////////////////////////////////////////////////
//
// Overlap function for a given x. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_x : public ComplexFunction
{
  public:

    Overlap_x(const SectionMode* m1_,const SectionMode* m2_,const Complex& x_)
      : m1(m1_), m2(m2_), x(x_) {}

    Complex operator()(const Complex& y)
    {
      counter++;

      Field f1 = m1->field(Coord(x,y,0));
      Field f2 = m2->field(Coord(x,y,0));

      return f1.E1 * f2.H2 - f1.E2 * f2.E1;
    }

  protected:

    const SectionMode* m1;
    const SectionMode* m2;
    Complex x;
};



/////////////////////////////////////////////////////////////////////////////
//
// Overlap function integrated out over y. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap : public ComplexFunction
{
  public:

    Overlap(const SectionMode* m1_, const SectionMode* m2_) 
      : m1(m1_), m2(m2_) {}

    Complex operator()(const Complex& x)
    {
      counter++;

      Overlap_x f(m1, m2, x);

      Wrap_real_to_real f_r(f);
      Wrap_real_to_imag f_i(f);

      // TODO: PML

      Real y_stop = real(m1->get_geom()->get_height());
      
      return patterson_quad(f_r, 0, y_stop, 1e-2, 4)
         + I*patterson_quad(f_i, 0, y_stop, 1e-2, 4);
    }

  protected:

    const SectionMode* m1;
    const SectionMode* m2;
};



/////////////////////////////////////////////////////////////////////////////
//
// overlap
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const SectionMode* mode_I,
                const SectionMode* mode_II,
                const SectionCache* cache,
                const std::vector<Complex>* disc,
                int i, int j, int I_index, int II_index)
{
  Overlap f(mode_I, mode_II);

  Wrap_real_to_real f_r(f);
  Wrap_real_to_imag f_i(f);

  // TODO: PML.

  Real x_stop = real(mode_I->get_geom()->get_width());

  return patterson_quad(f_r, 0, x_stop, 1e-2, 4)
     + I*patterson_quad(f_i, 0, x_stop, 1e-2, 4);
}
