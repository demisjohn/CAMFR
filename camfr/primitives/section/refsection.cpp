
/////////////////////////////////////////////////////////////////////////////
//
// File:     refsection.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20021104
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <algorithm>
#include "refsection.h"

using std::vector;
using std::sort;

#include "../../math/calculus/quadrature/patterson_quad.h"

/////////////////////////////////////////////////////////////////////////////
//
// Overlap function for a given x. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_x_ : public ComplexFunction
{
  public:

    Overlap_x_(const RefSectionMode* m1_,
               const RefSectionMode* m2_,
               const Section2D* profile_,
               const Complex& x_)
      : m1(m1_), m2(m2_), profile(profile_), x(x_) {}

    Complex operator()(const Complex& y)
    {
      counter++;

      Field f1 = m1->field(Coord(x,y,0));
      Field f2 = m2->field(Coord(x,y,0));

      Complex eps = profile->eps_at(Coord(x,y,0));

      //return f1.E1 * f2.H2 - f1.E2 * f2.H1; // Normalisation
      
      //return eps / eps0 * f1.Ez * f2.Ez; // O_zz

      return eps / eps0 * (f1.E1*f2.E1 + f1.E2*f2.E2); // O_trans
    }

  protected:

    const RefSectionMode* m1;
    const RefSectionMode* m2;
    const Section2D* profile;
    Complex x;
};



/////////////////////////////////////////////////////////////////////////////
//
// Overlap function integrated out over y. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_ : public ComplexFunction
{
  public:

    Overlap_(const RefSectionMode* m1_, const RefSectionMode* m2_,
             const Section2D* profile_, Slab* s_) 
      : m1(m1_), m2(m2_), profile(profile_), s(s_) {}

    Complex operator()(const Complex& x)
    {
      counter++;

      Overlap_x_ f(m1, m2, profile, x);

      Wrap_real_to_real f_r(f);
      Wrap_real_to_imag f_i(f);

      vector<Complex> slab_disc = s->get_discontinuities();

      // Loop over y materials.

      Complex result = 0.0;
      
      for (unsigned int l=0; l<slab_disc.size(); l++)
      {
        Complex y0 = l==0 ? 0.0 : slab_disc[l-1];
        Complex y1 = slab_disc[l];
      
        result += patterson_quad(f_r, real(y0), real(y1), 1e-2, 4)
          + I*patterson_quad(f_i, real(y0), real(y1), 1e-2, 4);
      }

      return result;
    }

  protected:

    const RefSectionMode* m1;
    const RefSectionMode* m2;
    const Section2D* profile;
    Slab* s;
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

Complex overlap_numeric_(const RefSectionMode* mode_I,
                         const RefSectionMode* mode_II,
                         const Section2D* profile)
{
  vector<Complex> section_disc = profile->discontinuities;

  Complex numeric = 0.0;

  for (unsigned int k=0; k<profile->slabs.size(); k++)
  {
    Complex x0 = k==0 ? 0.0 : section_disc[k-1];
    Complex x1 = section_disc[k];

    Slab* s = profile->slabs[k];

    Overlap_ f(mode_I, mode_II, profile, s);

    Wrap_real_to_real f_r(f);
    Wrap_real_to_imag f_i(f);

    numeric += patterson_quad(f_r, real(x0), real(x1), 1e-2, 4)
      + I*patterson_quad(f_i, real(x0), real(x1), 1e-2, 4);
  }
  
  return numeric;
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::find_modes
//  
/////////////////////////////////////////////////////////////////////////////

// TODO: check PML

void RefSection::find_modes()
{
  // Check values.

  if (global.lambda == 0)
  {
    py_error("Error: wavelength not set.");
    return;
  }
  
  if (global.N == 0)
  {
    py_error("Error: number of modes not set.");
    return;
  }

  Polarisation old_pol = global.polarisation;

  const int n = int(global.N/2);
  if (2*n != global.N)
    py_print("Warning: changing N to even number.");
  global.N = n;

  const Complex k2 = pow(2*pi/global.lambda * m->n(), 2);

  // Create TE and TM modes.
  
  Complex kx0 = (global_section.leftwall == E_wall) ? 0.0 : pi/2.;
  Complex ky0 = (global_slab.lowerwall   == NULL)   ? 0.0 : pi/2.;

  Complex x_offset = (global_section.rightwall == E_wall) ? 0.0 : pi/2.;
  Complex y_offset = (global_slab.upperwall    == NULL)   ? 0.0 : pi/2.;  

  vector<RefSectionMode*> TE_modes, TM_modes;
  
  for (int ix=0; ix<n; ix++)
    for (int iy=0; iy<n; iy++)
    {
      // Calc beta.

      Complex kx = (x_offset - kx0 + ix*pi)/a;
      Complex ky = (y_offset - ky0 + iy*pi)/b;

      if ((real(kx) < -1e-8) || (real(ky) < -1e-8))
        continue;

      Complex kt2 = kx*kx + ky*ky;
      Complex kz  = sqrt(k2 - kt2);

      if (real(kz) < 0)
        kz = -kz;

      if (abs(real(kz)) < 1e-12)
        if (imag(kz) > 0)
          kz = -kz;

      // Add TE mode?

      bool add_TE = true;
      
      if ( (abs(kx) < 1e-6) && (abs(ky) < 1e-6))
        add_TE = false;

      if ( (abs(kx) < 1e-6) && (abs(ky) > 1e-6))
        if (abs(cos(kx0)) < 1e-6)
          add_TE = false;

      if ( (abs(kx) > 1e-6) && (abs(ky) < 1e-6))
        if (abs(cos(ky0)) < 1e-6)
          add_TE = false;

      // Add TM mode?

      bool add_TM = true;
      
      if ( (abs(kx) < 1e-6) && (abs(ky) < 1e-6))
        add_TM = false;      

      if ( (abs(kx) < 1e-6) && (abs(ky) > 1e-6))
        if (abs(sin(kx0)) < 1e-6)
          add_TM = false;

      if ( (abs(kx) > 1e-6) && (abs(ky) < 1e-6))
        if (abs(sin(ky0)) < 1e-6)
          add_TM = false;

      // Add modes.

      if (add_TE)
      {
        RefSectionMode* mode = new RefSectionMode(TE,kz,kx,kx0,ky,ky0,this);
        TE_modes.push_back(mode);
      }

      if (add_TM)
      {
        RefSectionMode* mode = new RefSectionMode(TM,kz,kx,kx0,ky,ky0,this);
        TM_modes.push_back(mode);
      }
    }

  // Sort modes and create modeset.

  sort(TE_modes.begin(), TE_modes.end(), modesorter());
  sort(TM_modes.begin(), TM_modes.end(), modesorter());
  
  for (unsigned int i=n; i<TE_modes.size(); i++)
    delete TE_modes[i];  
  for (unsigned int i=n; i<TM_modes.size(); i++)
    delete TM_modes[i];

  TE_modes.erase(TE_modes.begin()+n, TE_modes.end());
  TM_modes.erase(TM_modes.begin()+n, TM_modes.end());

  for (unsigned int i=0; i<TE_modes.size(); i++)
    modeset.push_back(TE_modes[i]);
  for (unsigned int i=0; i<TM_modes.size(); i++)
    modeset.push_back(TM_modes[i]);

  //for (unsigned int i=0; i<modeset.size(); i++)
  //  std::cout << i 
  //            << overlap_numeric_(dynamic_cast<RefSectionMode*>(modeset[i]),
  //                                dynamic_cast<RefSectionMode*>(modeset[i])) 
  //            << std::endl;

  // Restore globals.

  global.polarisation = old_pol;
  global.N = 2*n;
}



/////////////////////////////////////////////////////////////////////////////
//
// cos_cos
//
//  Int(cos(ax+b).cos(cx+d).dx, x=x1..x2)
//  
/////////////////////////////////////////////////////////////////////////////

Complex cos_cos(const Complex& a,  const Complex& b,
                const Complex& c,  const Complex& d,
                const Complex& x1, const Complex& x2)
{
  Complex T1;
  if (abs(a+c) > 1e-6)
    T1 = 0.5/(a+c) * (sin((a+c)*x2+(b+d)) - sin((a+c)*x1+(b+d)));
  else
    T1 = 0.5*(x2-x1)*cos(b+d);

  Complex T2;
  if (abs(a-c) > 1e-6)
    T2 = 0.5/(a-c) * (sin((a-c)*x2+(b-d)) - sin((a-c)*x1+(b-d)));
  else
    T2 = 0.5*(x2-x1)*cos(b-d);

  return T1+T2;
}



/////////////////////////////////////////////////////////////////////////////
//
// sin_sin
//
//  Int(sin(ax+b).sin(cx+d).dx, x=x1..x2)
//  
/////////////////////////////////////////////////////////////////////////////

Complex sin_sin(const Complex& a,  const Complex& b,
                const Complex& c,  const Complex& d,
                const Complex& x1, const Complex& x2)
{
  Complex T1;
  if (abs(a+c) > 1e-6)
    T1 = 0.5/(a+c) * (sin((a+c)*x2+(b+d)) - sin((a+c)*x1+(b+d)));
  else
    T1 = 0.5*(x2-x1)*cos(b+d);

  Complex T2;
  if (abs(a-c) > 1e-6)
    T2 = 0.5/(a-c) * (sin((a-c)*x2+(b-d)) - sin((a-c)*x1+(b-d)));
  else
    T2 = 0.5*(x2-x1)*cos(b-d);

  return -T1+T2;
}



/////////////////////////////////////////////////////////////////////////////
 //
// RefSection::calc_overlap_matrices
//  
/////////////////////////////////////////////////////////////////////////////

void RefSection::calc_overlap_matrices(Section2D* profile,
                                       cMatrix* O_EE, cMatrix* O_MM,
                                       cMatrix* O_EM, cMatrix* O_zz)
{
  find_modes();

  const int n = int(global.N/2);

  vector<Complex> section_disc = profile->discontinuities;

  for (int i=1; i<=n; i++)
  {

    RefSectionMode* TE_i = dynamic_cast<RefSectionMode*>(modeset[i-1]);
    RefSectionMode* TM_i = dynamic_cast<RefSectionMode*>(modeset[n+i-1]);
    
    for (int j=i; j<=n; j++)
    {

      RefSectionMode* TE_j = dynamic_cast<RefSectionMode*>(modeset[j-1]);
      RefSectionMode* TM_j = dynamic_cast<RefSectionMode*>(modeset[n+j-1]);

      Complex O_EE_ij = 0.0; Complex O_MM_ij = 0.0;
      Complex O_EM_ij = 0.0; Complex O_ME_ij = 0.0;
      Complex O_zz_ij = 0.0;
     
      // Loop over x materials.
  
      for (unsigned int k=0; k<profile->slabs.size(); k++)
      {

        Complex x0 = k==0 ? 0.0 : section_disc[k-1];
        Complex x1 = section_disc[k];

        Slab* s = profile->slabs[k];
        vector<Complex> slab_disc = s->get_discontinuities();
    
        // Loop over y materials.

        for (unsigned int l=0; l<slab_disc.size(); l++)
        {
          Complex y0 = l==0 ? 0.0 : slab_disc[l-1];
          Complex y1 = slab_disc[l];

          Complex eps = s->eps_at(Coord(slab_disc[l],0,0,Min,Min,Min));

          // Calculate O_EE_ij.

          O_EE_ij += eps * TE_i->A * TE_j->A * TE_i->ky * TE_j->ky
            * cos_cos(TE_i->kx, TE_i->kx0, TE_j->kx, TE_j->kx0, x0, x1)
            * sin_sin(TE_i->ky, TE_i->ky0, TE_j->ky, TE_j->ky0, y0, y1)  +
            
            eps * TE_i->A * TE_j->A * TE_i->kx * TE_j->kx
            * sin_sin(TE_i->kx, TE_i->kx0, TE_j->kx, TE_j->kx0, x0, x1)
            * cos_cos(TE_i->ky, TE_i->ky0, TE_j->ky, TE_j->ky0, y0, y1);

          // Calculate O_MM_ij.

          O_MM_ij += eps * TM_i->A * TM_j->A * TM_i->kx * TM_j->kx
            * cos_cos(TM_i->kx, TM_i->kx0, TM_j->kx, TM_j->kx0, x0, x1)
            * sin_sin(TM_i->ky, TM_i->ky0, TM_j->ky, TM_j->ky0, y0, y1)  +
            
            eps * TM_i->A * TM_j->A * TM_i->ky * TM_j->ky
            * sin_sin(TM_i->kx, TM_i->kx0, TM_j->kx, TM_j->kx0, x0, x1)
            * cos_cos(TM_i->ky, TM_i->ky0, TM_j->ky, TM_j->ky0, y0, y1);

          // Calculate O_EM_ij.

          O_EM_ij += eps * TE_i->A * TM_j->A * TE_i->ky * TM_j->kx
            * cos_cos(TE_i->kx, TE_i->kx0, TM_j->kx, TM_j->kx0, x0, x1)
            * sin_sin(TE_i->ky, TE_i->ky0, TM_j->ky, TM_j->ky0, y0, y1)  -
            
            eps * TE_i->A * TM_j->A * TE_i->kx * TM_j->ky
            * sin_sin(TE_i->kx, TE_i->kx0, TM_j->kx, TM_j->kx0, x0, x1)
            * cos_cos(TE_i->ky, TE_i->ky0, TM_j->ky, TM_j->ky0, y0, y1);

          // Calculate O_ME_ij.

          O_ME_ij += eps * TM_i->A * TE_j->A * TM_i->kx * TE_j->ky
            * cos_cos(TM_i->kx, TM_i->kx0, TE_j->kx, TE_j->kx0, x0, x1)
            * sin_sin(TM_i->ky, TM_i->ky0, TE_j->ky, TE_j->ky0, y0, y1)  -
            
            eps * TM_i->A * TE_j->A * TM_i->ky * TE_j->kx
            * sin_sin(TM_i->kx, TM_i->kx0, TE_j->kx, TE_j->kx0, x0, x1)
            * cos_cos(TM_i->ky, TM_i->ky0, TE_j->ky, TE_j->ky0, y0, y1);

          // Calc O_zz_ij.

          O_zz_ij += -eps * TM_i->A * TM_j->A
            * TM_i->kt2() * TM_j->kt2() / TM_i->kz / TM_j->kz
            * sin_sin(TM_i->kx, TM_i->kx0, TM_j->kx, TM_j->kx0, x0, x1)
            * sin_sin(TM_i->ky, TM_i->ky0, TM_j->ky, TM_j->ky0, y0, y1);
        }
      }

      //Complex num = overlap_numeric_(TM_i,TE_j,profile);
      //std::cout << i << " " << j << " " << O_ME_ij/eps0 << " " 
      //          << num << std::endl;

      (*O_EE)(i,j) = O_EE_ij; (*O_EE)(j,i) = O_EE_ij;
      (*O_MM)(i,j) = O_MM_ij; (*O_MM)(j,i) = O_MM_ij;
      (*O_EM)(i,j) = O_EM_ij; (*O_EM)(j,i) = O_ME_ij; 
      (*O_zz)(i,j) = O_zz_ij; (*O_zz)(j,i) = O_zz_ij;
      
    }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::RefSectionMode
//  
/////////////////////////////////////////////////////////////////////////////

RefSectionMode::RefSectionMode(Polarisation pol,   const Complex& kz,
                               const Complex& kx_, const Complex& kx0_,
                               const Complex& ky_, const Complex& ky0_,
                               RefSection* geom_)
  : Mode(pol, kz, -kz), kx(kx_), kx0(kx0_), ky(ky_), ky0(ky0_), geom(geom_)
{
  // Normalise.

  const Complex kt2 = kx*kx + ky*ky;
  const Real omega  = 2*pi/global.lambda * c;

  const Complex C_TE = kz / omega / geom->get_mu();
  const Complex C_TM = omega * geom->get_eps() / kz;

  const Complex area = geom->get_width() * geom->get_height();

  Complex norm = (pol == TE) ? 0.25 * C_TE * kt2 * area
                             : 0.25 * C_TM * kt2 * area;

  if (abs(kx) < 1e-6)
    norm *= 2.0;
  
  if (abs(ky) < 1e-6)
    norm *= 2.0;

  if ((abs(kz) < 1e-6) || (abs(norm) < 1e-6))
    py_print("WARNING: reference mode close to cutoff!");

  A = 1.0/sqrt(norm);
}



/////////////////////////////////////////////////////////////////////////////
//
// RefSection::field
//  
/////////////////////////////////////////////////////////////////////////////

Field RefSectionMode::field(const Coord& coord) const
{
  // Calculate constants.

  const Complex kt2 = kx*kx + ky*ky;
  const Real omega  = 2*pi/global.lambda * c;

  const Complex C_TE = kz / omega / geom->get_mu();
  const Complex C_TM = omega * geom->get_eps() / kz;
 
  const Complex cos_x = cos(kx*coord.c1 + kx0); 
  const Complex sin_x = sin(kx*coord.c1 + kx0);

  const Complex cos_y = cos(ky*coord.c2 + ky0); 
  const Complex sin_y = sin(ky*coord.c2 + ky0);

  // Calculate field.

  Field f;

  f.E1 = (pol == TE) ?             ky*cos_x*sin_y :        kx*cos_x*sin_y;
  f.E2 = (pol == TE) ?            -kx*sin_x*cos_y :        ky*sin_x*cos_y;
  f.Ez = (pol == TE) ?                    0.      :  I*kt2/kz*sin_x*sin_y; 

  f.H1 = (pol == TE) ?        C_TE*kx*sin_x*cos_y :  -C_TM*ky*sin_x*cos_y;
  f.H2 = (pol == TE) ?        C_TE*ky*cos_x*sin_y :   C_TM*kx*cos_x*sin_y;
  f.Hz = (pol == TE) ? -I*kt2*C_TE/kz*cos_x*cos_y :               0.0;

  return A*f;
}

