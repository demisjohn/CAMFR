
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabmode.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000203
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "../../planar/planar.h"
#include "slaboverlap.h"
#include "slabwall.h"
#include "slabmode.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// SlabMode::get_kz()
//
/////////////////////////////////////////////////////////////////////////////

Complex SlabMode::get_kz() const 
{
  if (abs(global_slab.beta) < 1e-10)
    return kz;
  else // Rotate kz for off-axis incidence.
  {
    Complex new_kz = sqrt(kz*kz - pow(global_slab.beta, 2));

    if (imag(new_kz) > 0)
      new_kz = -new_kz;
    
    if (abs(imag(new_kz)) < abs(real(new_kz)))
      if (real(new_kz) < 0)
        new_kz = -new_kz;

    return new_kz;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabMode::field
//
/////////////////////////////////////////////////////////////////////////////

Field SlabMode::field(const Coord& coord_) const
{
  // Check and coerce input.

  Coord coord(coord_);
  coord.c1 = coord_.c1 + I*global_slab.left_PML;

  const Real x = real(coord.c1);
  const Real d = real(geom->get_width());

  if ( (x < 0) || (x > d) )
  {
    std::ostringstream s;
    s << "Error: x-value " << coord.c1 << " out of range.";
    py_error(s.str());
    exit (-1);
  }

  if (abs(x) < 1e-10)
    coord.c1_limit = Plus;

  if (abs(x-d) < 1e-10)
    coord.c1_limit = Min;
  
  // Calculate constants.
  
  const Complex k0 = 2*pi/global.lambda;
  const Complex kx = kx_at(coord);

  // Calculate amplitude of forward and backward waves at x-value.

  Complex right_x, left_x;
  forw_backw_at(coord, &right_x, &left_x);

  // Calculate total field.
  // Note that cs needs to be in sync with generalslab.cpp.
  
  const Complex sn = global_slab.beta / kz;
  const Complex cs = slab_signedsqrt(1.0 - sn*sn);
  
  Field field;

  if (pol == TE)
  {
    const Complex C = 1.0 / (k0*c) / geom->mu_at(coord);
        
    field.E1 = 0.0;
    field.E2 = right_x + left_x;
    field.Ez = 0.0;
    
    field.H1 = C * (-right_x - left_x) * kz;
    field.H2 = 0.0;
    field.Hz = C * ( right_x - left_x) * kx;

    if (abs(global_slab.beta) > 1e-6)
    {
      field.Ez  = -field.E2 * sn / sqrt(cs);
      field.E2 *= cs / sqrt(cs);
 
      field.H2  = field.Hz * sn / sqrt(cs);
      field.Hz *= cs / sqrt(cs);
      field.H1 /= sqrt(cs);
    }
  }
  else 
  { 
    const Complex C = 1.0 / (k0*c) / geom->eps_at(coord);
 
    field.H1 = 0.0;
    field.H2 = right_x - left_x;
    field.Hz = 0.0;

    field.E1 = C * ( right_x - left_x) * kz;
    field.E2 = 0.0;
    field.Ez = C * (-right_x - left_x) * kx;

    if (abs(global_slab.beta) > 1e-6)
    {
      field.E2  = field.Ez * sn / sqrt(cs);
      field.Ez *= cs / sqrt(cs);
      field.E1 /= sqrt(cs);
 
      field.Hz  = -field.H2 * sn / sqrt(cs);
      field.H2 *= cs / sqrt(cs);
    }
  }

  return field;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M_Mode::Slab_M_Mode
//
/////////////////////////////////////////////////////////////////////////////

Slab_M_Mode::Slab_M_Mode(Polarisation pol,   const Complex& kz, 
                         const Complex& kt_, const Slab_M* geom,
                         bool calc_fw)
    : SlabMode(pol, kz, geom), kt(kt_)
{
  A = (pol == TE) ? 0.0 : 1.0;
  B = (pol == TE) ? 1.0 : 0.0;

  calc_left_right(calc_fw);
}



/////////////////////////////////////////////////////////////////////////////
//
// safe_mult
//
/////////////////////////////////////////////////////////////////////////////

Complex safe_mult_(const Complex& a, const Complex& b)
{  
  Real threshold = pow(global.unstable_exp_threshold, 2);

  if ( (abs(a) > .1) && (abs(b) > .1) )
  {
    if (    (abs(a) < global.unstable_exp_threshold)
         || (abs(b) < global.unstable_exp_threshold) )
      return 0.0;
    else
      return a*b;
  }
  else
  {
    if (    (abs(a) < threshold * abs(b))
         || (abs(b) < threshold * abs(a)) )
      return 0.0;
    else
      return a*b;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M_Mode::calc_left_right
//
//  'right' refers to waves propagating in +x. 
//   'left' refers to waves propagating in -x.
//
//   The amplitudes of these waves can be calculated starting from the
//   first material (calc_fw == true) or from the last material
//   (calc_fw = false). 
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M_Mode::calc_left_right(bool calc_fw)
{ 
  const Slab_M* slab = dynamic_cast<const Slab_M*>(geom);
  
  const Complex k0_2 = pow(2*pi/global.lambda, 2);
  
  // Determine walls.

  SlabWall* l_wall = slab->leftwall  ? slab->leftwall  : global_slab.leftwall;
  SlabWall* r_wall = slab->rightwall ? slab->rightwall : global_slab.rightwall;
  
  // Calculate kx in each section.

  kx.clear();
  for (unsigned int i=0; i<slab->materials.size(); i++)
  {
    Material* mat = slab->materials[i];
    Complex n_2 = mat->epsr() * mat->mur();
    
    Complex kx_i = sqrt(k0_2*n_2 - kz*kz);

    if (abs(k0_2*n_2 - kz*kz) < 1e-10) // Improve stability for TEM mode.
      kx_i = 0.0;

    if (real(kx_i) < 0)
      kx_i = -kx_i;

    if (abs(real(kx_i)) < 1e-12)
      if (imag(kx_i) > 0)
        kx_i = -kx_i;
    
    kx.push_back(kx_i);
  }

  Complex R_wall_l = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_wall_r = r_wall ? r_wall->get_R12() : -1.0;

  Complex R_wall_inc = calc_fw ? R_wall_l : R_wall_r;
  Complex R_wall_ext = calc_fw ? R_wall_r : R_wall_l;
  
  // Calculate forward and backward plane wave expansion coefficients
  // at the interfaces.
  
  // Set seed field. If the left wall is not explicitly set, default to
  // electric wall. This also avoids a virtual function call.

  Complex fw_chunk_begin_scaled;
  Complex bw_chunk_begin_scaled;

  SlabWall* inc_wall = calc_fw ? l_wall : r_wall;

  if (!inc_wall)
  {
    fw_chunk_begin_scaled =  1.0;
    bw_chunk_begin_scaled = -1.0;
  }
  else
    inc_wall->get_start_field(&bw_chunk_begin_scaled, &fw_chunk_begin_scaled);

  if (calc_fw == false)
  {
    Complex swap = -fw_chunk_begin_scaled;
    fw_chunk_begin_scaled = -bw_chunk_begin_scaled;
    bw_chunk_begin_scaled = swap;
  }

  vector<Complex> fw, bw;
  
  fw.push_back(fw_chunk_begin_scaled);
  bw.push_back(bw_chunk_begin_scaled);
   
  // Loop through chunks and relate fields at the end of each chunk to
  // those at the beginning of the chunk.

  Complex fw_chunk_end_scaled;
  Complex bw_chunk_end_scaled;

  Complex scaling = 0.0;
  
  int k_start = calc_fw ?             0             : slab->materials.size()-1;
  int k_stop  = calc_fw ? slab->materials.size()-1  :             0;
  int delta_k = calc_fw ?             1             :            -1;
  
  for (int k=k_start; calc_fw ? k<=k_stop : k>=k_stop; k+=delta_k)
  {
    int i1 = (k==k_start) ? k_start : k-delta_k; // Index incidence medium.
    int i2 = k;                                  // Index exit medium.
    
    Material* mat1 = slab->materials[i1];
    Material* mat2 = slab->materials[i2];

    const Complex kx12 = ( (abs(kx[i1]) < 1e-10) && (abs(kx[i1]) < 1e-10) )
      ? 1.0 : kx[i1] / kx[i2];
    
    const Complex a = (pol == TE) ? kx12 * mat2->mu()  / mat1->mu()
                                  : kx12 * mat2->eps() / mat1->eps();
    
    const Real sign = (pol == TE) ? 1 : -1;
  
    // Cross the interface.
    
    fw_chunk_end_scaled =        (1.0+a)*0.5 * fw_chunk_begin_scaled +
                          sign * (1.0-a)*0.5 * bw_chunk_begin_scaled;

    bw_chunk_end_scaled = sign * (1.0-a)*0.5 * fw_chunk_begin_scaled +
                                 (1.0+a)*0.5 * bw_chunk_begin_scaled;

    // After a core region, reset the scaling.

    if ( real(mat1->eps_mu()) > real(mat2->eps_mu()) )
    {
      fw_chunk_end_scaled = safe_mult_(fw_chunk_end_scaled, exp(scaling));
      bw_chunk_end_scaled = safe_mult_(bw_chunk_end_scaled, exp(scaling));

      scaling = 0.0;
    }

    fw.push_back(safe_mult_(fw_chunk_end_scaled, exp(scaling)));
    bw.push_back(safe_mult_(bw_chunk_end_scaled, exp(scaling)));

    // Propagate in medium and scale along the way, by
    // factoring out and discarding the positive exponentials.
    // When R=0, we deal with an infinite medium.
    
    Complex I_kx_d = I * kx[i2] * slab->thicknesses[i2];

    if ( (k == k_start) && (abs(R_wall_inc) < 1e-10) ) 
      I_kx_d = 0;

    if ( (k == k_stop ) && (abs(R_wall_ext) < 1e-10) )
      I_kx_d = 0;

    if (real(I_kx_d) > 0.0)
    {
      fw_chunk_end_scaled *= exp(-2.0*I_kx_d);
      scaling += I_kx_d;
    }
    else
    {
      bw_chunk_end_scaled *= exp(+2.0*I_kx_d);
      scaling -= I_kx_d;
    }

    fw.push_back(safe_mult_(fw_chunk_end_scaled, exp(scaling)));
    bw.push_back(safe_mult_(bw_chunk_end_scaled, exp(scaling)));
    
    // Update values for next iteration.
    
    fw_chunk_begin_scaled = fw_chunk_end_scaled;
    bw_chunk_begin_scaled = bw_chunk_end_scaled;
  }

  // For metallic boundaries, correct for any rounding errors.

  if (real(kz) > imag(kz))
  {
    if (abs(R_wall_ext + 1.0) < 1e-10)
      bw.back() = -fw.back();
    if (abs(R_wall_ext - 1.0) < 1e-10)
      bw.back() =  fw.back();
  }

  // Assemble fw/bw into left/right. Note that the first element in
  // left/right is always duplicated.

   left.clear();
  right.clear();

  if (calc_fw)
  {
    right = fw;
     left = bw;
  }
  else
  {
    right.push_back(bw.back());
     left.push_back(fw.back());

    for (unsigned int i=bw.size()-1; i>0; i--)
    {
      right.push_back(bw[i]);
       left.push_back(fw[i]);
    } 
  }

  // For open boundary conditions, set the fields at infinity for the
  // leaky modes equal to zero to due to the infinite absorbion in the PML.

  if (abs(R_wall_l) < 1e-10)
  {
     left[0] =  left[1] = 0.0;
    right[0] = right[1] = 0.0;
  }
  
  if (abs(R_wall_r) < 1e-10)
  {
     left.back() = 0.0;
    right.back() = 0.0;
  }
  
  //for (unsigned int i=0; i<left.size(); i++)
  //  std::cout << i << kx[i] << left[i] << right[i] << std::endl;
  //std::cout << std::endl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M_Mode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M_Mode::normalise()
{
  Complex power = overlap(this, this);
  if (abs(power) < 1e-10)
  {
    py_print("Warning: mode close to cutoff.");
    power = 1;
  }
  
  A /= sqrt(power);
  B /= sqrt(power);

  for (unsigned int i=0; i<right.size(); i++)
  {
    right[i] /= sqrt(power);
     left[i] /= sqrt(power);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M_Mode::forw_backw_at
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M_Mode::forw_backw_at
  (const Coord& coord, Complex* fw, Complex* bw) const
{   
  // Calculate constants.
  
  const Complex kx = kx_at(coord);

  // Calculate index in chunk vector (first chunk is zero).
  
  const unsigned int index =
    index_lookup(coord.c1, coord.c1_limit, geom->discontinuities);

  // Calculate (positive) distance from previous discontinuity.
  
  const Complex d_prev =
    (index==0) ? coord.c1
               : coord.c1 - (geom->discontinuities)[index-1];
  
  // Calculate (positive) distance to next discontinuity.
  
  const Complex d_next = (geom->discontinuities)[index] - coord.c1;

  // Calculate indices in field vectors.

  const unsigned int prev_index = 2*index+1;
  const unsigned int next_index = 2*index+2;

  // If on a discontinuity, return without propagating.
  
  if (abs(d_prev) < 1e-9)
  {
    *fw = right[prev_index];
    *bw =  left[prev_index];

    return;
  }

  if (abs(d_next) < 1e-9)
  {
    *fw = right[next_index];
    *bw =  left[next_index];

    return;
  }
  
  // Else propagate only with decreasing exponentials.
  
  if (imag(kx) < 0) 
  {
    *fw = right[prev_index] * exp(-I * kx * d_prev);
    *bw =  left[next_index] * exp(-I * kx * d_next);
  }
  else
  {
    *fw = right[next_index] * exp( I * kx * d_next);
    *bw =  left[prev_index] * exp( I * kx * d_prev);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlabMode::UniformSlabMode
//
/////////////////////////////////////////////////////////////////////////////

UniformSlabMode::UniformSlabMode(Polarisation pol,
                                 const Complex& kz, const UniformSlab* geom)
  : SlabMode(pol, kz, geom)
{
  // Set A and B so that calc_interface gives results normalised in the same
  // way als not-uniform slabs.
  
  A = (pol == TE) ? 0.0 : 1./geom->get_core()->epsr();
  B = (pol == TE) ? 1./geom->get_core()->mur() : 0.0;

  // Set amplitude.

  SlabWall* l_wall = geom->leftwall ? geom->leftwall : global_slab.leftwall;
  
  if (!l_wall)
  {
    fw0 =  1;
    bw0 = -1;
  }
  else
    l_wall->get_start_field(&bw0, &fw0);
  
  // Set kx.
    
  const Complex k = 2*pi/global.lambda * geom->get_core()->n();
  
  kx = sqrt(k*k - kz*kz);

  if (abs(k-kz) < 1e-10) // Improve stability for TEM mode.
    kx = 0.0;

  if (real(kx) < 0)
    kx = -kx;

  if (abs(real(kx)) < 1e-12)
    if (imag(kx) > 0)
      kx = -kx;
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlabMode::normalise
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlabMode::normalise()
{
  Complex power = overlap(this, this);

  if (abs(power) < 1e-10)
  {
    py_print("Warning: mode close to cutoff.");
    power = 1;
  }

  A /= sqrt(power);
  B /= sqrt(power);
  
  fw0 /= sqrt(power);
  bw0 /= sqrt(power);
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlabMode::forw_backw_at
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlabMode::forw_backw_at
  (const Coord& coord, Complex* fw, Complex* bw) const
{ 
  // Propagate from x=0 forward to desired location.

  *fw = fw0 * exp(-I * kx * coord.c1);
  *bw = bw0 * exp( I * kx * coord.c1);
}
