
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

#include "../../../util/vectorutil.h"
#include "slaboverlap.h"

/////////////////////////////////////////////////////////////////////////////
//
// overlap
//  
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const SlabMode* mode_I,
                const SlabMode* mode_II,
                const SlabCache* cache=NULL,
                const vector<Complex>* disc=NULL,
                int i=0, int j=0,  int I_index=0, int II_index=0)
{
  // Check arguments.
  
  if (cache)
  {
    if ( (I_index < 1) || (I_index > 2) )
    {
      cerr << "Error: I_index should be 1 or 2." << endl;
      exit (-1);
    }
    if ( (II_index < 1) || (II_index > 2) )
    {
      cerr << "Error: II_index should be 1 or 2." << endl;
      exit (-1);
    }
  }
  
  // Set variables.

  const Real eps = 1e-10; // Don't choose too low.

  const SlabImpl* medium_I  = mode_I ->get_geom();
  const SlabImpl* medium_II = mode_II->get_geom();

  const Real omega = 2*pi/global.lambda * c;

  // TE and TM are orthogonal.

  if ((mode_I->pol) != (mode_II->pol))
    return 0.0;
  
  // Widths equal?
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    cerr << "Warning: complex widths don't match: "
         << medium_I ->get_width() << " and "
         << medium_II->get_width() << endl;
    exit (-1);
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

  Complex term1 = 0;
  Complex term2 = 0;

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
		   const SlabCache* cache=NULL,
		   const vector<Complex>* disc=NULL,
		   int i=0, int j=0, int I_index=0, int II_index=0)
{
  // Check arguments.
  
  if (cache)
  {
    if ( (I_index < 1) || (I_index > 2) )
    {
      cerr << "Error: I_index should be 1 or 2." << endl;
      exit (-1);
    }
    if ( (II_index < 1) || (II_index > 2) )
    {
      cerr << "Error: II_index should be 1 or 2." << endl;
      exit (-1);
    }
  }

  if ( (mode_I->pol != TM) || (mode_II->pol != TE) )
  {
     cerr << "Wrong polaristion of modes in overlap_TM_TE." << endl;
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
    cerr << "Warning: complex widths don't match: "
         << medium_I ->get_width() << " and "
         << medium_II->get_width() << endl;
    exit (-1);
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
