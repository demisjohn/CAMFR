
/////////////////////////////////////////////////////////////////////////////
//
// File:     slaboverlap.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000210
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "slaboverlap.h"
#include "slabmode.h"
#include "../../../math/calculus/quadrature/patterson_quad.h"

using std::vector;

#include "../../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// same
//
//  This function is needed to discrimate between calculating a true
//  normalisation integral and an integral between degenerate modes
//  in different cores.
//  
/////////////////////////////////////////////////////////////////////////////

bool same(const Complex& fw_I_l, const Complex& fw_II_l, 
          const Complex& bw_I_l, const Complex& bw_II_l,
          const Complex& fw_I_u, const Complex& fw_II_u, 
          const Complex& bw_I_u, const Complex& bw_II_u)
{
  const Real eps = 1e-4;

  if (abs(fw_I_l - fw_II_l) > eps)
    return false;
  if (abs(bw_I_l - bw_II_l) > eps)
    return false;
  if (abs(fw_I_u - fw_II_u) > eps)
    return false;
  if (abs(bw_I_u - bw_II_u) > eps)
    return false;

  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap
//  
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const SlabMode* mode_I,
                const SlabMode* mode_II,
                const SlabCache* cache,
                const vector<Complex>* disc,
                int i, int j, int I_index, int II_index)
{
  // Check arguments.
  
  if (cache)
  {
    if ( (I_index < 1) || (I_index > 2) )
    {
      py_error("Error: I_index should be 1 or 2.");
      exit (-1);
    }
    if ( (II_index < 1) || (II_index > 2) )
    {
      py_error("Error: II_index should be 1 or 2.");
      exit (-1);
    }
  }
  
  // Set variables.

  const Real eps = 1e-6; // Don't choose too low.

  const SlabImpl* medium_I  = mode_I ->get_geom();
  const SlabImpl* medium_II = mode_II->get_geom();

  const Real omega = 2*pi/global.lambda * c;

  // TE and TM are orthogonal.

  if ((mode_I->pol) != (mode_II->pol))
    return 0.0;
  
  // Widths equal?
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << medium_I ->get_width() << " and " << medium_II->get_width();
    py_error(s.str());
    return 0.0;
  }
  
  // Make sorted list of evaluation points for integrals (discontinuities).

  vector<Complex> local_disc;
  
  if (!disc)
  {
    local_disc = medium_I->discontinuities;

    local_disc.push_back(0.0);

    for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
      local_disc.push_back(medium_II->discontinuities[k]);

    remove_copies(&local_disc, 1e-9);

    sort(local_disc.begin(), local_disc.end(), RealSorter());

    disc=&local_disc;
  }
  
  // Actual calculation.

  Complex term1 = 0.0;
  Complex term2 = 0.0;  

  for (int k=0; k<int(disc->size()-1); k++) // Loop over all regions.
  {
    // Calculate necessary quantities in this region.

    const Complex d = (*disc)[k+1] - (*disc)[k];
    
    const Coord lower((*disc)[k],   0, 0, Plus);
    const Coord upper((*disc)[k+1], 0, 0, Min);

    const Complex kx_I  = mode_I ->kx_at(lower);
    const Complex kx_II = mode_II->kx_at(lower);

    const Complex C = (mode_I->pol == TE) ? 1.0 / medium_II->mu_at(lower)
                                          : 1.0 / medium_I->eps_at(lower);

    Complex fw_I_l, bw_I_l, fw_II_l, bw_II_l;
    Complex fw_I_u, bw_I_u, fw_II_u, bw_II_u;

    if (cache)
    {
       fw_I_l  = cache->fw_l(I_index,i,k+1);
       bw_I_l  = cache->bw_l(I_index,i,k+1);
       fw_I_u  = cache->fw_u(I_index,i,k+1);
       bw_I_u  = cache->bw_u(I_index,i,k+1);
   
       fw_II_l = cache->fw_l(II_index,j,k+1);
       bw_II_l = cache->bw_l(II_index,j,k+1);
       fw_II_u = cache->fw_u(II_index,j,k+1);
       bw_II_u = cache->bw_u(II_index,j,k+1);
    }
    else
    {
      mode_I ->forw_backw_at(lower, &fw_I_l,  &bw_I_l);    
      mode_I ->forw_backw_at(upper, &fw_I_u,  &bw_I_u);
      mode_II->forw_backw_at(lower, &fw_II_l, &bw_II_l);  
      mode_II->forw_backw_at(upper, &fw_II_u, &bw_II_u);
    }
    
    //
    // term 1 : Int(fw_I.fw_II + bw_I.bw_II)
    //
    
    if ( abs(kx_I + kx_II) > eps ) // normal case
    {
      // + f(upper_Min)

      term1 += C*I / (kx_I + kx_II) * ( fw_I_u * fw_II_u - bw_I_u * bw_II_u );
      
      // - f(lower_Plus)

      term1 -= C*I / (kx_I + kx_II) * ( fw_I_l * fw_II_l - bw_I_l * bw_II_l );
    }  
    else // normalisation integral (same modes) or degenerate case
    {
      // + f(upper_Min) - f(lower_Plus)

      if (same(fw_I_l, fw_II_l, bw_I_l, bw_II_l,
               fw_I_u, fw_II_u, bw_I_u, bw_II_u))
        term1 += C * d * ( fw_I_l * fw_II_l + bw_I_l * bw_II_l );
    }
    
    //
    // term 2 : Int(fw_I.bw_II + bw_I.fw_II)
    //
    
    if ( abs(kx_I - kx_II) > eps ) // normal case
    {
      // + f(upper_Min)
      
      term2 += C*I / (kx_I - kx_II) * ( fw_I_u * bw_II_u - bw_I_u * fw_II_u );
      
      // - f(lower_Plus)

      term2 -= C*I / (kx_I - kx_II) * ( fw_I_l * bw_II_l - bw_I_l * fw_II_l );
    }  
    else // normalisation integral (same modes) or degenerate case
    {
      // + f(upper_Min) - f(lower_Plus)

      if (same(fw_I_l, fw_II_l, bw_I_l, bw_II_l,
               fw_I_u, fw_II_u, bw_I_u, bw_II_u))
        term2 += C * d * ( fw_I_l * bw_II_l + bw_I_l * fw_II_l );
    }
        
  } // End loop over all regions.

  // Return final result.

  return (mode_I->pol == TE) ? mode_II->get_kz0() / omega * (term1 + term2)
                             : mode_I ->get_kz0() / omega * (term1 - term2);
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap_TM_TE
//  
/////////////////////////////////////////////////////////////////////////////

void overlap_TM_TE(const SlabMode* mode_I, const SlabMode* mode_II,
		   Complex* Ex_Hz, Complex* Ez_Hx,
		   const SlabCache* cache,
		   const vector<Complex>* disc,
		   int i, int j, int I_index, int II_index)
{
  // Check arguments.
  
  if (cache)
  {
    if ( (I_index < 1) || (I_index > 2) )
    {
      py_error("Error: I_index should be 1 or 2.");
      exit (-1);
    }
    if ( (II_index < 1) || (II_index > 2) )
    {
      py_error("Error: II_index should be 1 or 2.");
      exit (-1);
    }
  }

  if ( (mode_I->pol != TM) || (mode_II->pol != TE) )
  {
    py_error("Wrong polarisation of modes in overlap_TM_TE.");
    exit (-1);
  }

  // Set variables.

  const Real eps = 1e-10; // Don't choose too low.

  const SlabImpl* medium_I  = mode_I ->get_geom();
  const SlabImpl* medium_II = mode_II->get_geom();

  const Real omega = 2*pi/global.lambda * c;
  
  // Widths equal?
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << medium_I ->get_width() << " and " << medium_II->get_width();
    py_error(s.str());
    return;
  }
  
  // Make sorted list of evaluation points for integrals (discontinuities).

  vector<Complex> local_disc;
  
  if (!disc)
  {
    local_disc = medium_I->discontinuities;

    local_disc.push_back(0.0);

    for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
      local_disc.push_back(medium_II->discontinuities[k]);

    remove_copies(&local_disc, 1e-6);

    sort(local_disc.begin(), local_disc.end(), RealSorter());

    disc=&local_disc;
  }
  
  // Actual calculation.

  Complex term1_xz = 0.0; Complex term1_zx = 0.0;
  Complex term2_xz = 0.0; Complex term2_zx = 0.0;

  for (int k=0; k<int(disc->size()-1); k++) // Loop over all regions.
  {
    // Calculate necessary quantities in this region.

    const Complex d = (*disc)[k+1] - (*disc)[k];
    
    const Coord lower((*disc)[k],   0, 0, Plus);
    const Coord upper((*disc)[k+1], 0, 0, Min);

    const Complex kx_I  = mode_I ->kx_at(lower);
    const Complex kx_II = mode_II->kx_at(lower);

    Complex C_xz = kx_II / (medium_I->eps_at(lower) * medium_II->mu_at(lower));
    Complex C_zx = kx_I  / (medium_I->eps_at(lower) * medium_II->mu_at(lower));

    Complex fw_I_l, bw_I_l, fw_II_l, bw_II_l;
    Complex fw_I_u, bw_I_u, fw_II_u, bw_II_u;

    if (cache)
    {
       fw_I_l  = cache->fw_l(I_index,i,k+1);
       bw_I_l  = cache->bw_l(I_index,i,k+1);
       fw_I_u  = cache->fw_u(I_index,i,k+1);
       bw_I_u  = cache->bw_u(I_index,i,k+1);
   
       fw_II_l = cache->fw_l(II_index,j,k+1);
       bw_II_l = cache->bw_l(II_index,j,k+1);
       fw_II_u = cache->fw_u(II_index,j,k+1);
       bw_II_u = cache->bw_u(II_index,j,k+1);
    }
    else
    {
      mode_I ->forw_backw_at(lower, &fw_I_l,  &bw_I_l);    
      mode_I ->forw_backw_at(upper, &fw_I_u,  &bw_I_u);
      mode_II->forw_backw_at(lower, &fw_II_l, &bw_II_l);  
      mode_II->forw_backw_at(upper, &fw_II_u, &bw_II_u);
    }
    
    //
    // term 1 : Int(fw_I.fw_II + bw_I.bw_II)
    //

    Complex inc1 = 0.0;
    
    if ( abs(kx_I + kx_II) > eps ) // normal case
    {
      // + f(upper_Min)
      
      inc1 += I / (kx_I + kx_II) * ( fw_I_u * fw_II_u - bw_I_u * bw_II_u );

      // - f(lower_Plus)

      inc1 -= I / (kx_I + kx_II) * ( fw_I_l * fw_II_l - bw_I_l * bw_II_l );
    }  
    else // normalisation integral (same modes) or degenerate case
    {
      // + f(upper_Min) - f(lower_Plus)

      inc1 += d * ( fw_I_l * fw_II_l + bw_I_l * bw_II_l );
    }

    term1_xz += C_xz * inc1; 
    term1_zx += C_zx * inc1;
    
    //
    // term 2 : Int(fw_I.bw_II + bw_I.fw_II)
    //
    
    Complex inc2 = 0.0;
    
    if ( abs(kx_I - kx_II) > eps ) // normal case
    {
      // + f(upper_Min)
      
      inc2 += I / (kx_I - kx_II) * ( fw_I_u * bw_II_u - bw_I_u * fw_II_u );
      
      // - f(lower_Plus)

      inc2 -= I / (kx_I - kx_II) * ( fw_I_l * bw_II_l - bw_I_l * fw_II_l );
    }  
    else // normalisation integral (same modes) or degenerate case
    {
      // + f(upper_Min) - f(lower_Plus)

      inc2 += d * ( fw_I_l * bw_II_l + bw_I_l * fw_II_l );
    }

    term2_xz += C_xz * inc2; 
    term2_zx += C_zx * inc2;
        
  } // End loop over all regions.

  // Return final result.

  *Ex_Hz =  mode_I->get_kz0() / omega / omega * (term1_xz - term2_xz);
  *Ez_Hx = mode_II->get_kz0() / omega / omega * (term1_zx + term2_zx);
}



/////////////////////////////////////////////////////////////////////////////
//
// Overlap reference modes
//
//  Numeric integration, only used for verification purposes. 
//
/////////////////////////////////////////////////////////////////////////////

class Overlap_ : public ComplexFunction
{
  public:

    Overlap_(const SlabMode* m1_,const SlabMode* m2_,
             const Slab_M* profile_)
      : m1(m1_), m2(m2_), profile(profile_){}

    Complex operator()(const Complex& x)
    {
      counter++;

      Field f1 = m1->field(Coord(x,0,0));
      Field f2 = m2->field(Coord(x,0,0));

      Complex eps = profile->eps_at(Coord(x,0,0));
      
      return eps / eps0 * f1.Ez * f2.Ez; // O_zz

      //return eps / eps0 * (f1.E1*f2.E1 + f1.E2*f2.E2); // O_trans
    }

  protected:

    const SlabMode* m1;
    const SlabMode* m2;
    const Slab_M*   profile;
};



/////////////////////////////////////////////////////////////////////////////
//
// overlap_numeric
//
//   Verification function which calculates the overlap integrals 
//   numerically. Slower and less acurate.
//   To be removed after the analytical overlap calculation matures.
//
/////////////////////////////////////////////////////////////////////////////

Complex overlap_numeric_(const SlabMode* mode_I,
                         const SlabMode* mode_II,
                         const Slab_M* profile)
{
  vector<Complex> slab_disc = profile->get_discontinuities();

  Complex numeric = 0.0;

  for (unsigned int k=0; k<slab_disc.size(); k++)
  {
    Complex x0 = k==0 ? 0.0 : slab_disc[k-1];
    Complex x1 = slab_disc[k];

    Overlap_ f(mode_I, mode_II, profile);

    Wrap_real_to_real f_r(f);
    Wrap_real_to_imag f_i(f);

    numeric += patterson_quad(f_r, real(x0), real(x1), 1e-6, 4)
      + I*patterson_quad(f_i, real(x0), real(x1), 1e-6, 4);
  }
  
  return numeric;
}



/////////////////////////////////////////////////////////////////////////////
//
// overlap_reference_modes
//
/////////////////////////////////////////////////////////////////////////////

void overlap_reference_modes(cMatrix* O_tt, cMatrix* O_zz,
                             const UniformSlab& ref, const Slab_M& profile)
{
  // Set constants.

  const Real eps_copies = 1e-10; // Don't choose too low.
  const Complex C_TM = global.lambda / (2.*pi*c * ref.eps_at(Coord(0,0,0)));
  
  // Widths equal?
  
  if (abs(ref.get_width() - profile.get_width()) > eps_copies)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << ref.get_width() << " and " << profile.get_width();
    py_error(s.str());
    return;
  }
  
  // Make list of evaluation points for integrals (discontinuities).

  vector<Complex> disc = profile.get_discontinuities();
  disc.insert(disc.begin(), 0.0);
  
  // Actual calculation.

  for (int k=0; k<int(disc.size()-1); k++) // Loop over all regions.
  {
    const Complex d = disc[k+1] - disc[k];

    const Coord lower(disc[k],   0, 0, Plus);
    const Coord upper(disc[k+1], 0, 0, Min);

    const Complex eps = profile.eps_at(lower);

    for (int i=1; i<=ref.N(); i++)
    {
      UniformSlabMode* mode_I 
        = dynamic_cast<UniformSlabMode*>(ref.get_mode(i));

      Complex fw_I_l, bw_I_l, fw_I_u, bw_I_u;

      mode_I->forw_backw_at(lower, &fw_I_l,  &bw_I_l);
      mode_I->forw_backw_at(upper, &fw_I_u,  &bw_I_u);

      const Complex kz_I = mode_I->get_kz0();
      const Complex kx_I = mode_I->kx_at(lower);
    
      for (int j=i; j<=ref.N(); j++)
      {
        UniformSlabMode* mode_II
          = dynamic_cast<UniformSlabMode*>(ref.get_mode(j));

        Complex fw_II_l, bw_II_l, fw_II_u, bw_II_u;

        mode_II->forw_backw_at(lower, &fw_II_l, &bw_II_l);
        mode_II->forw_backw_at(upper, &fw_II_u, &bw_II_u);

        const Complex kz_II = mode_II->get_kz0();
        const Complex kx_II = mode_II->kx_at(lower);
    
        //
        // term 1 : Int(fw_I.fw_II + bw_I.bw_II)
        //
    
        Complex term1 = 0.0;
        if ( abs(kx_I + kx_II) > eps_copies ) // normal case
        {
          // + f(upper_Min)

          term1 += I/(kx_I + kx_II) * ( fw_I_u * fw_II_u - bw_I_u * bw_II_u );
      
          // - f(lower_Plus)

          term1 -= I/(kx_I + kx_II) * ( fw_I_l * fw_II_l - bw_I_l * bw_II_l );
        }  
        else // normalisation integral (same modes) or degenerate case
        {
          // + f(upper_Min) - f(lower_Plus)

          term1 += d * ( fw_I_l * fw_II_l + bw_I_l * bw_II_l );      
        }
    
        //
        // term 2 : Int(fw_I.bw_II + bw_I.fw_II)
        //
  
        Complex term2 = 0.0;  
        if ( abs(kx_I - kx_II) > eps_copies ) // normal case
        {
          // + f(upper_Min)
      
          term2 += I/(kx_I - kx_II) * ( fw_I_u * bw_II_u - bw_I_u * fw_II_u );
      
          // - f(lower_Plus)

          term2 -= I/(kx_I - kx_II) * ( fw_I_l * bw_II_l - bw_I_l * fw_II_l );
        }  
        else // normalisation integral (same modes) or degenerate case
        {
          // + f(upper_Min) - f(lower_Plus)

          term2 += d * ( fw_I_l * bw_II_l + bw_I_l * fw_II_l );
        }

        //
        // Update matrices.
        //

        if (global.polarisation == TE)
          (*O_tt)(i,j) += eps * (term1 + term2);
        else
        {
          (*O_tt)(i,j) += eps * C_TM*C_TM * (term1 - term2) * kz_I*kz_II;
          (*O_zz)(i,j) += eps * C_TM*C_TM * (term1 + term2) * kx_I*kx_II;
        }
      }
    }
  }

  // Fill out symmetric part of the matrices.

  for (int i=1; i<=ref.N(); i++)
    for (int j=1; j<i; j++)
    {
      (*O_tt)(i,j) = (*O_tt)(j,i);

      if (global.polarisation == TM)
        (*O_zz)(i,j) = (*O_zz)(j,i);      
    }

  // Verify numerically.

/*
  for (int i=1; i<=ref.N(); i++)
  {
    UniformSlabMode* mode_I
      = dynamic_cast<UniformSlabMode*>(ref.get_mode(i));

    for (int j=1; j<=i; j++)
    {
      UniformSlabMode* mode_II
        = dynamic_cast<UniformSlabMode*>(ref.get_mode(j));

      Complex num = overlap_numeric_(mode_I,mode_II,&profile);
      if (abs(num - (*O_zz)(i,j)/eps0) > 1e-3)
      {
        std::cout << abs(mode_I->get_kz() - mode_II->get_kz()) << " " ;
        std::cout << i << " " << j << " " << (*O_zz)(i,j)/eps0 << " "
                  << num << std::endl;
      }    
    } 
  }
  
  exit(-1);
*/
}

