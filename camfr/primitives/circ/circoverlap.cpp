
/////////////////////////////////////////////////////////////////////////////
//
// File:     circoverlap.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990519
// Version:  1.0
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <algorithm>
#include "../../util/vectorutil.h"
#include "circoverlap.h"
#include "circmode.h"

/////////////////////////////////////////////////////////////////////////////
//
// overlap
//  
/////////////////////////////////////////////////////////////////////////////

Complex overlap(const Circ_M_Mode* mode_I,
                const Circ_M_Mode* mode_II,
                const CircCache*   cache,
                const vector<Complex>* disc,
                int i, int j, int I_index, int II_index)
{ 
  // Check arguments.

  if (cache)
  {
    if ( (I_index < 1) || (I_index > 2) )
    {
      cerr << "Error: I_index should be 1 or 2." << endl;
      exit(-1);
    }
    if ( (II_index < 1) || (II_index > 2) )
    {
      cerr << "Error: II_index should be 1 or 2." << endl;
      exit(-1);
    }
  }

  // Set variables.

  const Real eps = 1e-10; // Don't choose too low.
  
  const Real n = global_circ.order;

  const Circ_M* medium_I  = mode_I ->get_geom();
  const Circ_M* medium_II = mode_II->get_geom();

  const Real omega = 2*pi/global.lambda * c;
  
  // Outer radii equal?
  
  if (medium_I->radius[medium_I->M-1] != medium_II->radius[medium_II->M-1])
  {
    cerr << "Error: outer radii don't match: "
         << medium_I ->radius[medium_I ->M-1] << " and "
         << medium_II->radius[medium_II->M-1] << endl;
    exit (-1);
  }
  
  // Make sorted list of evaluation points for integrals (discontinuities).

  vector<Complex> local_disc;

  if (!disc)
  {
    local_disc = medium_I->radius;

    local_disc.push_back(0.0);

    for (unsigned int k=0; k<medium_II->radius.size(); k++)
      local_disc.push_back(medium_II->radius[k]);

    remove_copies(&local_disc, 1e-6);

    sort(local_disc.begin(), local_disc.end(), RealSorter());

    disc = &local_disc;
  }
  
  // Actual calculation.

  Complex term1 = 0;
  Complex term2 = 0;
  Complex term3 = 0;
  Complex term4 = 0;

  for (int k=0; k<int(disc->size()-1); k++) // Loop over all regions.
  {   
    // Calculate necessary quantities in this region.

    const Complex r_l = (*disc)[k];   // lower
    const Complex r_u = (*disc)[k+1]; // upper
    
    const Coord lower(r_l,0,0, Plus);
    const Coord upper(r_u,0,0, Min);

    const Complex kr_I_2  = pow(mode_I ->kr_at(lower), 2);
    const Complex kr_II_2 = pow(mode_II->kr_at(lower), 2);

    const Complex mu_I   = medium_I -> mu_at(lower);
    const Complex eps_II = medium_II->eps_at(lower);

    Complex  E_I_l,  H_I_l,  E_II_l,  H_II_l;
    Complex dE_I_l, dH_I_l, dE_II_l, dH_II_l;
    Complex  E_I_u,  H_I_u,  E_II_u,  H_II_u;
    Complex dE_I_u, dH_I_u, dE_II_u, dH_II_u;

    if (cache)
    {
       E_I_l  = cache-> E_l(I_index,i,k+1);
       H_I_l  = cache-> H_l(I_index,i,k+1);
       E_I_u  = cache-> E_u(I_index,i,k+1);
       H_I_u  = cache-> H_u(I_index,i,k+1);
      dE_I_l  = cache->dE_l(I_index,i,k+1);
      dH_I_l  = cache->dH_l(I_index,i,k+1);
      dE_I_u  = cache->dE_u(I_index,i,k+1);
      dH_I_u  = cache->dH_u(I_index,i,k+1);

       E_II_l = cache-> E_l(II_index,j,k+1);
       H_II_l = cache-> H_l(II_index,j,k+1);
       E_II_u = cache-> E_u(II_index,j,k+1);
       H_II_u = cache-> H_u(II_index,j,k+1);
      dE_II_l = cache->dE_l(II_index,j,k+1);
      dH_II_l = cache->dH_l(II_index,j,k+1);
      dE_II_u = cache->dE_u(II_index,j,k+1);
      dH_II_u = cache->dH_u(II_index,j,k+1);
    }
    else
    {
      // Factor out the part of the field that has an angular dependence
      // with an integrated value of n*pi over 0..2*pi.

      const bool ang_dep = false;
    
      Field f_I_l  = mode_I ->field_at(lower, &dE_I_l,  &dH_I_l,  ang_dep);
      Field f_II_l = mode_II->field_at(lower, &dE_II_l, &dH_II_l, ang_dep);
      Field f_I_u  = mode_I ->field_at(upper, &dE_I_u,  &dH_I_u,  ang_dep);
      Field f_II_u = mode_II->field_at(upper, &dE_II_u, &dH_II_u, ang_dep);

      E_I_l  =  f_I_l.Ez; H_I_l  =  f_I_l.Hz;
      E_II_l = f_II_l.Ez; H_II_l = f_II_l.Hz;
      E_I_u  =  f_I_u.Ez; H_I_u  =  f_I_u.Hz;
      E_II_u = f_II_u.Ez; H_II_u = f_II_u.Hz;
    } 

    //
    // term 1 (..kzi.kzj.gradEi x gradHj)
    //
    
    if (n != 0) // speedup
    {
      const Complex C1 = 1.0 / kr_I_2 / kr_II_2;

      // + f(upper_Min) - f(lower_Plus)

      term1 += C1 * (   E_I_u * H_II_u
                      - E_I_l * H_II_l );
    }
    
    //
    // term 2 (..kzi.epsj.gradEi . gradEj)
    //

    const Complex C2 = 1.0 / kr_I_2 / kr_II_2 * eps_II;
    
    if ( abs(kr_I_2 - kr_II_2) > eps ) // normal case
    {
      // + f(upper_Min)
      
      term2 += C2 * r_u / (kr_I_2 - kr_II_2) * (   kr_I_2  * E_I_u  * dE_II_u
                                                 - kr_II_2 * E_II_u * dE_I_u);

      // - f(lower_Plus)
      
      term2 -= C2 * r_l / (kr_I_2 - kr_II_2) * (   kr_I_2  * E_I_l  * dE_II_l
                                                 - kr_II_2 * E_II_l * dE_I_l);
    } 
    else // normalisation integral (same modes) or degenerate case
    {
      // + f(upper_Min)
      
      term2 += C2 * (                 r_u            *  E_I_u * dE_II_u
                      + 0.5 * (r_u*r_u*kr_I_2 - n*n) *  E_I_u *  E_II_u
                      + 0.5 *        r_u*r_u         * dE_I_u * dE_II_u );

      // - f(lower_Plus)
      
      term2 -= C2 * (                 r_l            *  E_I_l * dE_II_l
                      + 0.5 * (r_l*r_l*kr_I_2 - n*n) *  E_I_l *  E_II_l
                      + 0.5 *        r_l*r_l         * dE_I_l * dE_II_l );
    }
    
    //
    // term 3 (..kzj.mui.gradHi . gradHj)
    //

    const Complex C3 = 1.0 / kr_I_2 / kr_II_2 * mu_I;
    
    if ( abs(kr_I_2 - kr_II_2) > eps ) // normal case
    {
      
      // + f(upper_Min)
      
      term3 += C3 * r_u / (kr_I_2 - kr_II_2) * (   kr_I_2  * H_I_u  * dH_II_u
                                                 - kr_II_2 * H_II_u * dH_I_u);

      // - f(lower_Plus)
      
      term3 -= C3 * r_l / (kr_I_2 - kr_II_2) * (   kr_I_2  * H_I_l  * dH_II_l
                                                 - kr_II_2 * H_II_l * dH_I_l);
    }
    else // normalisation integral (same modes) or degenerate case
    { 
      // + f(upper_Min)
      
      term3 += C3 * (                 r_u            *  H_I_u * dH_II_u
                      + 0.5 * (r_u*r_u*kr_I_2 - n*n) *  H_I_u *  H_II_u
                      + 0.5 *        r_u*r_u         * dH_I_u * dH_II_u );

      // - f(lower_Plus)
      
      term3 -= C3 * (                 r_l            *  H_I_l * dH_II_l
                      + 0.5 * (r_l*r_l*kr_I_2 - n*n) *  H_I_l *  H_II_l
                      + 0.5 *        r_l*r_l         * dH_I_l * dH_II_l );
    }
    
    //
    // term 4 (..mui.epsj.gradHi x gradEj)
    //
    
    if (n != 0) // speedup
    {
      const Complex C4 = 1.0 / kr_I_2 / kr_II_2 * mu_I * eps_II;

      // + f(upper_Min) - f(lower_Plus)

      term4 += C4 * (   H_I_u * E_II_u
                      - H_I_l * E_II_l );
    }
    
  } // end loop over all regions
  
  // Return final result.

  const Real ang_integral = (n == 0) ? 2*pi : pi;
  const Real sign = (global_circ.fieldtype == cos_type) ? 1 : -1;

  const Complex kz_I  = mode_I ->get_kz();
  const Complex kz_II = mode_II->get_kz();  
  
  return -ang_integral * (   term1 * kz_I  * kz_II * n *  sign
                           + term2 * kz_I  * omega
                           + term3 * kz_II * omega
                           - term4 * omega * omega * n * -sign );
}
