
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabdisp.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
//
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "../../planar/planar.h"
#include "slabwall.h"
#include "slabdisp.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// SlabDisp::SlabDisp
//
/////////////////////////////////////////////////////////////////////////////

SlabDisp::SlabDisp(const Expression& expression, Real lambda_,
                   SlabWall* lowerwall_, SlabWall* upperwall_)
  : lambda(lambda_), lowerwall(lowerwall_), upperwall(upperwall_)
{
  // Create a table with the eps's and the mu's. This avoids repeated
  // virtual function calls and allows to sweep these parameters without
  // tampering with the original materials.

  material_expression_to_table(expression, &eps, &mu, &thicknesses);
    
  // Find min refractive index of structure.

  min_eps_mu = eps[0]*mu[0];
  
  for (unsigned int i=1; i<eps.size(); i++)
  {
    Complex eps_mu = eps[i]*mu[i];

    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabDisp::SlabDisp
//
/////////////////////////////////////////////////////////////////////////////

SlabDisp::SlabDisp(const vector<Material*>& materials,
                   const vector<Complex>& thicknesses_, Real lambda_,
                   SlabWall* lowerwall_, SlabWall* upperwall_)
  : thicknesses(thicknesses_), lambda(lambda_),
    lowerwall(lowerwall_), upperwall(upperwall_)
{
  // Create a table with the eps's and the mu's. This avoids repeated
  // virtual function calls and allows to sweep these parameters without
  // tampering with the original materials.

  for (unsigned int i=0; i<materials.size(); i++)
  {
    eps.push_back(materials[i]->eps());
     mu.push_back(materials[i]->mu());
  }

  // Find min refractive index of structure.

  min_eps_mu = eps[0]*mu[0];
  
  for (unsigned int i=1; i<eps.size(); i++)
  {
    Complex eps_mu = eps[i]*mu[i];
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabDisp::operator()
//
//   kt = kx in the region with minimum refractive index.
//
/////////////////////////////////////////////////////////////////////////////

Complex SlabDisp::operator()(const Complex& kt)
{
  counter++;

  global.lambda = lambda;

  const Complex C = pow(2*pi/lambda, 2) / (eps0 * mu0);

  // Calculate kx in all materials.

  vector<Complex> kx;
  
  for (unsigned int i=0; i<eps.size(); i++)
  {
    Complex kx_i = sqrt(C*(eps[i]*mu[i] - min_eps_mu) + kt*kt);

    if (real(kx_i) < 0)
      kx_i = -kx_i;

    if (abs(real(kx_i)) < 1e-12)
      if (imag(kx_i) > 0)
        kx_i = -kx_i;

    kx.push_back(kx_i);
  }

  // For SlabWall_PC, we still need to set Planar::kt.

  Complex beta = sqrt(C*min_eps_mu - kt*kt);

  if (real(beta) < 0)
    beta = -beta;
  
  if (abs(imag(beta)) < 1e-12)
    if (imag(beta) > 0)
      beta = -beta;

  Planar::set_kt(beta);

  // Determine wall reflectivities. If the walls are not explicitly set,
  // default to electric walls. This also avoids a virtual function call.

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;
  
  // Set seed field.
  
  Complex fw_chunk_begin_scaled;
  Complex bw_chunk_begin_scaled;

  if (!l_wall)
  {
    fw_chunk_begin_scaled =  1.0;
    bw_chunk_begin_scaled = -1.0;
  }
  else
    l_wall->get_start_field(&bw_chunk_begin_scaled, &fw_chunk_begin_scaled);
  
  // Loop through chunks and relate fields at the end of each chunk to
  // those at the beginning of the chunk.
  // A chunk consists of crossing an interface and propagation in the
  // exit medium.

  Complex fw_chunk_end_scaled;
  Complex bw_chunk_end_scaled;
  
  for (unsigned int k=0; k<eps.size(); k++)
  {
    unsigned int i1 = (k==0) ? 0 : k-1; // Index incidence medium.
    unsigned int i2 = k;                // Index exit medium.
    
    const Complex a = (global.polarisation == TE)
      ? kx[i1] / kx[i2] *  mu[i2] /  mu[i1]
      : kx[i1] / kx[i2] * eps[i2] / eps[i1];
    
    const Real sign = (global.polarisation == TE) ? 1 : -1;
  
    // Cross the interface.

    fw_chunk_end_scaled =        (1.0+a)*0.5 * fw_chunk_begin_scaled +
                          sign * (1.0-a)*0.5 * bw_chunk_begin_scaled;

    bw_chunk_end_scaled = sign * (1.0-a)*0.5 * fw_chunk_begin_scaled +
                                 (1.0+a)*0.5 * bw_chunk_begin_scaled;

    // Propagate in medium and scale along the way, by
    // factoring out and discarding the positive exponentials.
    // When R=0, we deal with an infinite medium.
    
    Complex I_kx_d = I * kx[i2] * thicknesses[i2];

    if ( (k == 0) && (abs(R_lower) < 1e-10) ) 
      I_kx_d = 0;

    if ( (k == eps.size()-1) && (abs(R_upper) < 1e-10) )
      I_kx_d = 0;

    if (global.solver == track)
    {
      if (real(I_kx_d) > 0)
        fw_chunk_end_scaled *= exp(-2.0*I_kx_d);
      else
        bw_chunk_end_scaled *= exp(+2.0*I_kx_d);
    }
    else // No scaling to keep the function analytic.
    {
      fw_chunk_end_scaled *= exp(-I_kx_d);      
      bw_chunk_end_scaled *= exp(+I_kx_d);
    }
    
    // Update values for next iteration.
    
    fw_chunk_begin_scaled = fw_chunk_end_scaled;
    bw_chunk_begin_scaled = bw_chunk_end_scaled;
  }

  // Return error.
  
  if (!u_wall)
    return fw_chunk_end_scaled + bw_chunk_end_scaled;
  else
    return u_wall->get_error(fw_chunk_end_scaled, bw_chunk_end_scaled);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabDisp::get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> SlabDisp::get_params() const
{
  // Extract params from pointers to thicknesses and to materials in
  // expression.

  vector<Complex> params;
  
  for (unsigned int i=0; i<thicknesses.size(); i++)
  {
    params.push_back(thicknesses[i]);
    params.push_back(eps[i]);
    params.push_back(mu[i]);
  }

  params.push_back(lambda);

  params.push_back(min_eps_mu);

  return params;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabDisp::set_params
//
/////////////////////////////////////////////////////////////////////////////

void SlabDisp::set_params(const vector<Complex>& params)
{
  unsigned int params_index = 0;
    
  for (unsigned int i=0; i<thicknesses.size(); i++)
  {
    thicknesses[i] = params[params_index++];
    eps[i]         = params[params_index++];
    mu[i]          = params[params_index++];
  }

  lambda = real(params[params_index++]);

  min_eps_mu = params[params_index];
}
