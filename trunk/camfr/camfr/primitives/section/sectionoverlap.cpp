
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
#include "../slab/isoslab/slabmode.h"
#include "../slab/isoslab/slaboverlap.h"
#include "../../math/calculus/quadrature/patterson_quad.h"

using std::vector;

#include "../../util/vectorutil.h"

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

      return f1.E1 * f2.H2 - f1.E2 * f2.H1;
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
  // Numeric calculation.

  Overlap f(mode_I, mode_II);

  Wrap_real_to_real f_r(f);
  Wrap_real_to_imag f_i(f);

  // TODO: PML.

  Real x_stop = real(mode_I->get_geom()->get_width());

  Complex numeric = patterson_quad(f_r, 0, x_stop, 1e-2, 4)
     + I*patterson_quad(f_i, 0, x_stop, 1e-2, 4);

  return numeric;
}



/////////////////////////////////////////////////////////////////////////////
//
// signedsqrt2__
//
//   Square root with branch cut at 45 degrees.
//
/////////////////////////////////////////////////////////////////////////////

Complex signedsqrt2__(const Complex& kz2)
{
  Complex new_kz = sqrt(kz2);

  if (imag(new_kz) > 0)
    new_kz = -new_kz;

  if (abs(imag(new_kz)) < abs(real(new_kz)))
    if (real(new_kz) < 0)
      new_kz = -new_kz;

  return new_kz;
}



/////////////////////////////////////////////////////////////////////////////
//
// Int(exp(k.z), z=0..d)
//
/////////////////////////////////////////////////////////////////////////////

inline Complex int_exp(const Complex& k, const Complex& d)
{
  if (abs(k) < 1e-10)
    return d;
  else
    return (exp(k*d)-1.0)/k;
}



/////////////////////////////////////////////////////////////////////////////
//
// z_integral
//
/////////////////////////////////////////////////////////////////////////////

Complex z_integral
  (const Complex& fw_I,  const Complex& bw_I,  SlabMode* mode_I,
   const Complex& fw_II, const Complex& bw_II, SlabMode* mode_II,
   const Complex& kz_I,  const Complex& kz_II, const Complex& d, bool minus)
{
  // Determine betas.

  Complex old_beta = global_slab.beta;
  
  global_slab.beta = kz_I;
  Complex j_beta_I = I*mode_I->get_kz();
  
  global_slab.beta = kz_II;
  Complex j_beta_II = I*mode_II->get_kz();

  global_slab.beta = old_beta;

  Complex s = minus ? -1.0 : 1.0;

  // Calculate integral.

  return       fw_I * fw_II * int_exp(-j_beta_I - j_beta_II, d)
         + s * fw_I * bw_II * int_exp(-j_beta_I + j_beta_II, d)
         + s * bw_I * fw_II * int_exp( j_beta_I - j_beta_II, d)
         +     bw_I * bw_II * int_exp( j_beta_I + j_beta_II, d);
}




/////////////////////////////////////////////////////////////////////////////
//
// overlap_slice
//
/////////////////////////////////////////////////////////////////////////////

// TODO: Sectioncache for get_fw_bw.

void overlap_slice(const Section2D& sec_I, const Section2D& sec_II, 
                   const Slab& slab_I,     const Slab& slab_II,
                   const Complex& z_start, const Complex& z_stop,
                   cMatrix* O)
{
  std::cout << "overlap between" 
            << *slab_I.get_core() << " " << *slab_II.get_core()
            << z_start << z_stop << std::endl;
  
  const Complex d = z_stop-z_start;
  
  for (int i=1; i<=global.N; i++)
  {
    std::cout << "overlap " << i << std::endl << std::flush;
    
    cVector fw_I(sec_I.get_M(),fortranArray);
    cVector bw_I(sec_I.get_M(),fortranArray); 

    dynamic_cast<SectionMode*>(sec_I.get_mode(i))
      ->get_fw_bw(z_start, Plus, &fw_I, &bw_I);
    
    for (int j=1; j<=global.N; j++)
    {
      cVector fw_II(sec_II.get_M(),fortranArray);
      cVector bw_II(sec_II.get_M(),fortranArray); 

      dynamic_cast<SectionMode*>(sec_II.get_mode(j))
        ->get_fw_bw(z_start, Plus, &fw_II, &bw_II);

      //
      // Term 1: Ez_I * H1_II
      //

      Complex term1 = 0.0;
      for (int jj=1; jj<=sec_II.get_M(); jj++)
      {
        SlabMode* mode_II = dynamic_cast<SlabMode*>(slab_II.get_mode(jj));
        Complex kz_II = sec_II.get_mode(j)->get_kz();

        Complex sn_II = kz_II / mode_II->get_kz0();
        Complex cs_II = signedsqrt2__(1.0 - sn_II*sn_II);
        
        if (slab_II.get_mode(jj)->pol != TM)
          for (int ii=1; ii<=sec_I.get_M(); ii++)
          {
            SlabMode* mode_I = dynamic_cast<SlabMode*>(slab_I.get_mode(ii));
            Complex kz_I = sec_I.get_mode(i)->get_kz();

            Complex sn_I = kz_I / mode_I->get_kz0();
            Complex cs_I = signedsqrt2__(1.0 - sn_I*sn_I);

            // Integration over x.

            Complex x_factor;
            if (mode_I->pol == TE)
              x_factor = sn_I/sqrt(cs_II)/sqrt(cs_I)*overlap(mode_I,mode_II);
            else
            {
              Complex E1_Hz, Ez_H1;
              overlap_TM_TE(mode_I, mode_II, &E1_Hz, &Ez_H1);
              x_factor = cs_I/sqrt(cs_II)/sqrt(cs_I) * Ez_H1;
            }

            // Integration over z.

            Complex z_factor = z_integral(fw_I (ii),bw_I (ii),mode_I,
                                          fw_II(jj),bw_II(jj),mode_II, 
                                          kz_I,kz_II,d,true);

            term1 += x_factor*z_factor;
          }
      }
      (*O)(i,j) += term1;



      //
      // Term 2: -E1_I * Hz_II
      //

      Complex term2 = 0.0;
      for (int ii=1; ii<=sec_I.get_M(); ii++)
      {
        SlabMode* mode_I = dynamic_cast<SlabMode*>(slab_I.get_mode(ii));
        Complex kz_I = sec_I.get_mode(i)->get_kz();
        
        Complex sn_I = kz_I / mode_I->get_kz0();
        Complex cs_I = signedsqrt2__(1.0 - sn_I*sn_I);

        if (mode_I->pol != TE)
          for (int jj=1; jj<=sec_II.get_M(); jj++)
          {
            SlabMode* mode_II = dynamic_cast<SlabMode*>(slab_II.get_mode(jj));
            Complex kz_II = sec_II.get_mode(j)->get_kz();

            Complex sn_II = kz_II / mode_II->get_kz0();
            Complex cs_II = signedsqrt2__(1.0 - sn_II*sn_II);

            // Integration over x.

            Complex x_factor;
            if (mode_II->pol == TE)
            {
              Complex E1_Hz, Ez_H1;
              overlap_TM_TE(mode_I, mode_II, &E1_Hz, &Ez_H1);
              x_factor = cs_II/sqrt(cs_II)/sqrt(cs_I) * E1_Hz;
            }
            else
              x_factor = -sn_II/sqrt(cs_II)/sqrt(cs_I)*overlap(mode_I,mode_II);
            
            // Integration over z.

            Complex z_factor = z_integral(fw_I (ii),bw_I (ii),mode_I,
                                          fw_II(jj),bw_II(jj),mode_II, 
                                          kz_I,kz_II,d,false);

            term2 += x_factor*z_factor;
          }
      }

      (*O)(i,j) -= term2; 
    } 
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap_matrices
//
/////////////////////////////////////////////////////////////////////////////

void overlap_matrices
  (cMatrix* O, SectionImpl* medium_I_, SectionImpl* medium_II_)
{

  Section2D* medium_I  = dynamic_cast<Section2D*>(medium_I_);
  Section2D* medium_II = dynamic_cast<Section2D*>(medium_II_);
  
  // Widths equal?

  const Real eps = 1e-10; // Don't choose too low.
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << medium_I ->get_width() << " and " << medium_II->get_width();
    py_error(s.str());
    exit (-1);
  }
  
  // Make sorted list of separation points between different slabs.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);

  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  sort(disc.begin(), disc.end(), RealSorter());

  *O = 0.0;
  for (unsigned int k=0; k<disc.size()-1; k++)
    overlap_slice(*medium_I, *medium_II,
                  *medium_I ->slabs[index_lookup(disc[k], Plus,
                                                 medium_I ->discontinuities)],
                  *medium_II->slabs[index_lookup(disc[k], Plus,
                                                 medium_II->discontinuities)],
                  disc[k], disc[k+1], O);    
}

