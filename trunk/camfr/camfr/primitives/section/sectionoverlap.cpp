
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
#include "../slab/slabmatrixcache.h"
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

      // Make list of intervals.

      const SectionImpl* sec_I  
        = dynamic_cast<SectionImpl*>(m1->get_geom());
      const SectionImpl* sec_II 
        = dynamic_cast<SectionImpl*>(m2->get_geom());

      const Slab& slab_I  
        = *sec_I ->slabs[index_lookup(x, Plus, sec_I ->discontinuities)];

      const Slab& slab_II 
        = *sec_II->slabs[index_lookup(x, Plus, sec_II->discontinuities)];

      vector<Complex> disc = slab_I.get_discontinuities();

      disc.insert(disc.begin(), 0.0);

      vector<Complex> disc_II = slab_II.get_discontinuities();
      for (unsigned int k=0; k<disc_II.size(); k++)
        disc.push_back(disc_II[k]);

      remove_copies(&disc, 1e-9);

      sort(disc.begin(), disc.end(), RealSorter());

      // Do numeric integration.

      Overlap_x f(m1, m2, x);

      Wrap_real_to_real f_r(f);
      Wrap_real_to_imag f_i(f);

      Complex result = 0.0;
      for (unsigned int k=0; k<disc.size()-1; k++)
      {
        Real y_start = real(disc[k]);
        Real y_stop = real(disc[k+1]);

        result += patterson_quad(f_r, y_start, y_stop, 1e-2, 4)
              + I*patterson_quad(f_i, y_start, y_stop, 1e-2, 4);
      }
      
      return result;
    }

  protected:

    const SectionMode* m1;
    const SectionMode* m2;
};



/////////////////////////////////////////////////////////////////////////////
//
// overlap_numeric
//
//   Verification function which calculates the overlap integrals 
//   numerically. Slower and less acurate.
//   TODO: check PML in z direction.
//   To be removed after the analytical overlap calculation matures.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap_numeric(const SectionMode* mode_I,
                        const SectionMode* mode_II)
{
  // Make sorted list of separation points between different slabs.

  Section2D* medium_I  = dynamic_cast<Section2D*>(mode_I ->get_geom());
  Section2D* medium_II = dynamic_cast<Section2D*>(mode_II->get_geom());

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);

  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);
  
  sort(disc.begin(), disc.end(), RealSorter());

  // Loop over slices.

  Overlap f(mode_I, mode_II);

  Wrap_real_to_real f_r(f);
  Wrap_real_to_imag f_i(f);

  Complex result = 0.0;
  for (unsigned int k=0; k<disc.size()-1; k++)
  {
    Real x_start = real(disc[k]);
    Real x_stop  = real(disc[k+1]);

    result += patterson_quad(f_r, x_start, x_stop, 1e-2, 4)
          + I*patterson_quad(f_i, x_start, x_stop, 1e-2, 4);
  }

  return result;
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

  Complex old_beta = global.slab_ky;
  
  global.slab_ky = kz_I;
  Complex j_beta_I = I*mode_I->get_kz();
  
  global.slab_ky = kz_II;
  Complex j_beta_II = I*mode_II->get_kz();

  global.slab_ky = old_beta;

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

Complex overlap_slice(const SectionMode* sec_I_mode, 
                      const SectionMode* sec_II_mode,
                      const Complex& z_start, const Complex& z_stop,
                      FieldExpansion* field_I, FieldExpansion* field_II,
                      OverlapMatrices* m, int I_index, int II_index)
{
  Complex old_beta = global.slab_ky;

  // Calculate field expansion at start of slice.

  const Complex kz_I  = sec_I_mode ->get_kz();
  const Complex kz_II = sec_II_mode->get_kz();
  
  const Complex d = z_stop-z_start;

  const SectionImpl* sec_I  
    = dynamic_cast<SectionImpl*>( sec_I_mode->get_geom());
  const SectionImpl* sec_II 
    = dynamic_cast<SectionImpl*>(sec_II_mode->get_geom());

  const Slab& slab_I  
    = *sec_I ->slabs[index_lookup(z_start, Plus, sec_I ->discontinuities)];

  const Slab& slab_II 
    = *sec_II->slabs[index_lookup(z_start, Plus, sec_II->discontinuities)];

  const int sec_I_M  = sec_I ->get_M2();
  const int sec_II_M = sec_II->get_M2();

  // Get field expansion from cache or calculate it from scratch.

  cVector* fw_I; cVector* bw_I;
  if (field_I)
  {
    fw_I = &(field_I->fw);
    bw_I = &(field_I->bw);
  }
  else
  {
    fw_I = new cVector(sec_I_M, fortranArray);
    bw_I = new cVector(sec_I_M, fortranArray);

    sec_I_mode->get_fw_bw(z_start, Plus, fw_I, bw_I);
  } 

  cVector* fw_II; cVector* bw_II;
  if (field_II)
  {
    fw_II = &(field_II->fw);
    bw_II = &(field_II->bw);
  }
  else
  {
    fw_II = new cVector(sec_II_M, fortranArray);
    bw_II = new cVector(sec_II_M, fortranArray);

    sec_II_mode->get_fw_bw(z_start, Plus, fw_II, bw_II);
  }

  //
  // Term 1: Ez_I * H1_II
  //

  Complex term1 = 0.0;
  for (int jj=1; jj<=sec_II_M; jj++)
  { 
    SlabMode* mode_II = dynamic_cast<SlabMode*>(slab_II.get_mode(jj));

    global.slab_ky = sec_II_mode->get_kz();

    Complex sn_II = mode_II->get_sin();
    Complex cs_II = mode_II->get_cos();

    if (slab_II.get_mode(jj)->pol != TM)
      for (int ii=1; ii<=sec_I_M; ii++)
      {
        SlabMode* mode_I = dynamic_cast<SlabMode*>(slab_I.get_mode(ii));

        global.slab_ky = sec_I_mode->get_kz();

        Complex sn_I = mode_I->get_sin();
        Complex cs_I = mode_I->get_cos();

        // Integration over x.

        Complex x_factor;
        if (mode_I->pol == TE) // Same polarisation.
        {
          Complex o;
          if (mode_I == mode_II)
            o = (ii == jj) ? 1.0 : 0.0;
          else
          {
            if (!m)
              o = overlap(mode_I,mode_II);
            else
              o = m->TE_TE(I_index,ii,jj); 
          }
          x_factor = sn_I/sqrt(cs_II)/sqrt(cs_I) * o;
        }
        else // Cross polarisation.
        {
          if (!m)
          {
            Complex Ex_Hz, Ez_Hx;
            overlap_TM_TE(mode_I, mode_II, &Ex_Hz, &Ez_Hx);
            x_factor = cs_I/sqrt(cs_II)/sqrt(cs_I) * Ez_Hx;
          }
          else
            x_factor = cs_I/sqrt(cs_II)/sqrt(cs_I) 
                            * m->Ez_Hx_cross(I_index,ii-sec_I_M/2,jj);
        }

        // Integration over z.

        Complex z_factor = z_integral((*fw_I )(ii),(*bw_I )(ii),mode_I,
                                      (*fw_II)(jj),(*bw_II)(jj),mode_II,
                                      kz_I,kz_II,d,true);

        term1 += x_factor*z_factor;
      }
  }
      
  //
  // Term 2: -E1_I * Hz_II
  //

  Complex term2 = 0.0;
  for (int ii=1; ii<=sec_I_M; ii++)
  {
    SlabMode* mode_I = dynamic_cast<SlabMode*>(slab_I.get_mode(ii));
 
    global.slab_ky = sec_I_mode->get_kz();

    Complex sn_I = mode_I->get_sin();
    Complex cs_I = mode_I->get_cos();

    if (mode_I->pol != TE)
      for (int jj=1; jj<=sec_II_M; jj++)
      {
        SlabMode* mode_II = dynamic_cast<SlabMode*>(slab_II.get_mode(jj));
    
        global.slab_ky = sec_II_mode->get_kz();

        Complex sn_II = mode_II->get_sin();
        Complex cs_II = mode_II->get_cos();

        // Integration over x.

        Complex x_factor;
        if (mode_II->pol == TE) // Cross polarisation.
        {
          if (!m)
          {
            Complex Ex_Hz, Ez_Hx;
            overlap_TM_TE(mode_I, mode_II, &Ex_Hz, &Ez_Hx);
            x_factor = cs_II/sqrt(cs_II)/sqrt(cs_I) * Ex_Hz;
          }
          else
            x_factor = cs_II/sqrt(cs_II)/sqrt(cs_I) 
                            * m->Ex_Hz_cross(I_index,ii-sec_I_M/2,jj);
        }
        else // Same polarisation.
        {
          Complex o;         
          if (mode_I == mode_II)
            o = (ii == jj) ? 1.0 : 0.0;
          else
          {
            if (!m)
              o = overlap(mode_I,mode_II);
            else
              o = m->TM_TM(I_index,ii-sec_I_M/2,jj-sec_II_M/2);
          }
          x_factor = -sn_II/sqrt(cs_II)/sqrt(cs_I) * o;
        }
        
        // Integration over z.

        Complex z_factor = z_integral((*fw_I )(ii),(*bw_I )(ii),mode_I,
                                      (*fw_II)(jj),(*bw_II)(jj),mode_II,
                                      kz_I,kz_II,d,false);

        term2 += x_factor*z_factor;
      }
  }

  // Clean up and return result.

  if (!field_I)
  {
    delete fw_I;  delete bw_I;
    delete fw_II, delete bw_II;
  }

  global.slab_ky = old_beta;

  return term1 - term2;
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap_pw
//
//   Overlap between plane wave based modes.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap_pw(const Section2D_Mode* sec_I_mode, 
                   const Section2D_Mode* sec_II_mode)
{
  Complex result = 0.0;

  const int M = global_section.M;
  const int N = global_section.N;

  const Complex W = sec_I_mode->get_geom()->get_width();
  const Complex H = sec_I_mode->get_geom()->get_height();

  cVector* Ex = sec_I_mode ->Ex;
  cVector* Ey = sec_I_mode ->Ey;
  cVector* Hx = sec_II_mode->Hx;
  cVector* Hy = sec_II_mode->Hy;

  for (Real m=-M; m<=M; m+=1.0)
    for (Real n=-N; n<=N; n+=1.0)
    {
      int i1 = int(( m+M+1) + ( n+N)*(2*M+1));
      int i2 = int((-m+M+1) + (-n+N)*(2*M+1));

      result += (*Ex)(i1)*(*Hy)(i2) - (*Ey)(i1)*(*Hx)(i2);
    }

  //Complex numeric = overlap_numeric(sec_I_mode, sec_II_mode);

  //std::cout << "overlap " << numeric << " " << W*H*result
  //           << " " << numeric-W*H*result << std::endl;
  
  //if (abs(numeric - W*H*result) > 1e-4)
  //  std::cout << "****" << std::endl;
  
  return W*H*result;
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap
//
//   General overlap function.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const Section2D_Mode* sec_I_mode,
                const Section2D_Mode* sec_II_mode)
{
  Section2D* medium_I  = dynamic_cast<Section2D*>(sec_I_mode ->geom);
  Section2D* medium_II = dynamic_cast<Section2D*>(sec_II_mode->geom);

  //
  // Two uncorrected modes.
  //

  if (!sec_I_mode->corrected && !sec_II_mode->corrected)
    return overlap_pw(sec_I_mode, sec_II_mode);



  //
  // Two corrected modes.
  //
  
  if (sec_I_mode->corrected && sec_II_mode->corrected)
  {
    // Make sorted list of separation points between different slabs.

    vector<Complex> disc = medium_I->discontinuities;

    disc.push_back(0.0);

    for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
      disc.push_back(medium_II->discontinuities[k]);

    remove_copies(&disc, 1e-6);

    sort(disc.begin(), disc.end(), RealSorter());

    // Loop over slices.

    Complex result = 0.0;
    for (unsigned int k=0; k<disc.size()-1; k++)
      result += overlap_slice(sec_I_mode, sec_II_mode, disc[k], disc[k+1]);
    return result;
  }



  //
  // Corrected and uncorrected mode.
  //

  const Section2D* medium = (sec_I_mode->corrected) ? medium_I : medium_II;

  const Section2D_Mode* sec_mode = (sec_I_mode->corrected) ? sec_I_mode
                                                           : sec_II_mode;

  vector<Complex> disc = medium->discontinuities;
  disc.insert(disc.begin(), 0.0);

  const int M = global_section.M;
  const int N = global_section.N;

  int sec_M = sec_mode->get_geom()->get_M2();

  cVector fw(sec_M, fortranArray);
  cVector bw(sec_M, fortranArray);

  cVector* Ex = sec_I_mode ->Ex;
  cVector* Ey = sec_I_mode ->Ey;
  cVector* Hx = sec_II_mode->Hx;
  cVector* Hy = sec_II_mode->Hy;

  const Complex W = sec_mode->get_geom()->get_width();
  const Complex H = sec_mode->get_geom()->get_height();

  int old_N = global.N;
  global.N = sec_M;  

  int old_orthogonal = global.orthogonal;
  global.orthogonal = false;

  Complex old_beta = global.slab_ky;
  global.slab_ky = sec_mode->get_kz();
  
  Complex result = 0.0;

  for (unsigned int k=0; k<disc.size()-1; k++) // Loop over slices.
  {
    SlabImpl* slab = medium->
      slabs[index_lookup(disc[k],Plus,medium->discontinuities)]->get_impl();

    sec_mode->get_fw_bw(disc[k], Plus, &fw, &bw);

    const Complex x0 = disc[k];
    const Complex x1 = disc[k+1];

    const Complex d = x1 - x0;

    for (int jj=1; jj<=sec_M; jj++)
    { 
      const SlabMode* slab_mode = dynamic_cast<SlabMode*>(slab->get_mode(jj));
      const Complex kz = slab_mode->get_kz();

      for (Real m=-M; m<=M; m+=1.0)
      {
        // TODO: alpha0, beta0

        Complex alpha = 2.*pi*(m/W/2.);

        Complex fw_jj = fw(jj)*int_exp(-I*kz+I*alpha,d)*exp(I*alpha*x0);
        Complex bw_jj = bw(jj)*int_exp(+I*kz+I*alpha,d)*exp(I*alpha*x0);
        
        for (Real n=-N; n<=N; n+=1.0)
        {
          int i1 = int((m+M+1) + (n+N)*(2*M+1));
          Complex beta = 2.*pi*(n/H/2.);
          Complex Ox, Oy;

          // TODO: precalc overlap_pw.

          if (sec_I_mode->corrected)
          {
            overlap_pw(slab_mode, beta, true, &Ox, &Oy);
            result += (fw_jj-bw_jj)*Ox*(*Hy)(i1) - (fw_jj+bw_jj)*Oy*(*Hx)(i1);
          }
          else
          {            
            overlap_pw(slab_mode, beta, false, &Ox, &Oy);
            result += (*Ex)(i1)*(fw_jj-bw_jj)*Oy - (*Ey)(i1)*(fw_jj+bw_jj)*Ox;
          }
        }
      }
    }
  }

  global.N = old_N;
  global.slab_ky = old_beta;
  global.orthogonal = old_orthogonal;

  return result;
  
  // Numerical verification.

  Complex numeric = overlap_numeric(sec_I_mode, sec_II_mode);
  
  if (abs(numeric - result) > 1e-3)
  {
    std::cout << sec_I_mode->corrected << " " << sec_II_mode->corrected 
              << " " << sec_I_mode->n_eff() << " " << sec_II_mode->n_eff()
              << std::endl;
    std::cout << "overlap " << numeric << " " << result
              << " " << abs(numeric-result) <<  " " 
              << abs((result)/numeric) << std::endl;
  }

  return result;
}




