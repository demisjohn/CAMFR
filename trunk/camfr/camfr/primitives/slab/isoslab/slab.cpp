
/////////////////////////////////////////////////////////////////////////////
//
// File:     slab.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
//
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "slab.h"
#include "slabwall.h"
#include "slabdisp.h"
#include "slabmode.h"
#include "slaboverlap.h"
#include "../slabmatrixcache.h"
#include "../../../math/calculus/calculus.h"

using std::vector;

#include "../../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::Slab_M
//
/////////////////////////////////////////////////////////////////////////////

Slab_M::Slab_M(const Expression& expression)
{
  Expression ex = expression.flatten();
  Complex current_x = 0.0;
  
  for (unsigned int i=0; i<ex.get_size(); i++)
  {
    Material* m = dynamic_cast<Material*>(ex.get_term(i)->get_mat());

    if (!m)
    {
      py_error("Error: expression contains non-material term.");
      exit (-1);
    }

    Complex thickness = ex.get_term(i)->get_d();

    // Combine two succesive terms containing the same material.

    if ( (i+1 < ex.get_size()) && (m == ex.get_term(i+1)->get_mat()) )
    {
      thickness += ex.get_term(i+1)->get_d();
      i++;
    }

    if (abs(thickness))
    {
      current_x += thickness;
      
      materials.push_back(m);
      thicknesses.push_back(thickness);
      discontinuities.push_back(current_x);
    }
  }

  // Determine core.

  unsigned int core_index = 0;
  for (unsigned int i=0; i<materials.size(); i++)
    if ( real(materials[i]->eps_mu()) > real(materials[core_index]->eps_mu()) )
      core_index = i;

  core = materials[core_index];
  uniform = false;  
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::Slab_M
//
/////////////////////////////////////////////////////////////////////////////

Slab_M::Slab_M(const Slab_M& slab)
{
  params          = slab.params;
  materials       = slab.materials;
  thicknesses     = slab.thicknesses;
  discontinuities = slab.discontinuities;
  core            = slab.core;
  uniform         = slab.uniform;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::operator=
//
/////////////////////////////////////////////////////////////////////////////

Slab_M& Slab_M::operator=(const Slab_M& slab)
{
  if (this == &slab)
    return *this;

  params          = slab.params;
  materials       = slab.materials;
  thicknesses     = slab.thicknesses;
  discontinuities = slab.discontinuities;
  core            = slab.core;
  uniform         = slab.uniform;
  
  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::operator==
//
//   Contrary to other waveguides, we will be lazy here and compare the
//   adresses instead of the actual layouts...
//
/////////////////////////////////////////////////////////////////////////////

bool Slab_M::operator==(const Waveguide& w) const
{
  return (this == &w);
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool Slab_M::no_gain_present() const
{
  for (unsigned int i=0; i<materials.size(); i++)
    if (!materials[i]->no_gain_present())
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::eps_at
//
/////////////////////////////////////////////////////////////////////////////

Complex Slab_M::eps_at(const Coord& coord) const
{
  unsigned int i = index_lookup(coord.c1, coord.c1_limit, discontinuities);
  
  if (i == materials.size())
  {
    py_error("Warning: coordinate out of range in eps_at. Restricting it.");
    i = materials.size()-1;
  }

  return materials[i]->eps();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::mu_at
//
/////////////////////////////////////////////////////////////////////////////

Complex Slab_M::mu_at(const Coord& coord) const
{
  unsigned int i = index_lookup(coord.c1, coord.c1_limit, discontinuities);
  
  if (i == materials.size())
  {
    py_error("Warning: coordinate out of range in mu_at. Restricting it.");
    i = materials.size()-1;
  }

  return materials[i]->mu();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes()
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

  if (!recalc_needed())
    return;


  // Only TE or TM modes needed.

  if ((global.polarisation == TE) || (global.polarisation == TM))
  {
    vector<Complex> old_kt;
    for (unsigned int i=0; i<modeset.size(); i++)
      old_kt.push_back(dynamic_cast<Slab_M_Mode*>(modeset[i])->get_kt());

    vector<Complex> kt(find_kt(old_kt));

    build_modeset(kt);  
  }

  // Both TE and TM modes needed.

  if (global.polarisation == TE_TM)
  {
    // Cheat on global variables.

    const int n = int(global.N/2);
    if (2*n != global.N)
      py_print("Warning: changing N to even number.");
    global.N = n;

    const Complex old_beta = global_slab.beta;
    global_slab.beta = 0.0;

    vector<Complex> old_params = params;

    // Find TE modes.

    global.polarisation = TE;

    vector<Complex> old_kt_TE;
    if (modeset.size())
      for (unsigned int i=0; i<n; i++)
        old_kt_TE.push_back(dynamic_cast<Slab_M_Mode*>(modeset[i])->get_kt()); 

    vector<Complex> kt(find_kt(old_kt_TE));

    // Find TM modes

    params = old_params;

    last_lambda = 0.0; // Force a recalc.

    global.polarisation = TM;

    vector<Complex> old_kt_TM;
    if (modeset.size())
      for (unsigned int i=n; i<2*n; i++)
        old_kt_TM.push_back(dynamic_cast<Slab_M_Mode*>(modeset[i])->get_kt());

    vector<Complex> kt_TM(find_kt(old_kt_TM));

    // Restore global variables and build modeset.

    global.N = 2*n;
    global_slab.beta = old_beta;
    global.polarisation = TE_TM;

    kt.insert(kt.end(), kt_TM.begin(), kt_TM.end());
    build_modeset(kt);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::find_kt(vector<Complex>& old_kt)
{
  // If we already calculated modes for a different wavelength/gain
  // combination, use these as an initial estimate, else find them
  // from scratch.

  vector<Complex> kt;

  if (global.sweep_from_previous && (modeset.size() >= global.N))
      kt = find_kt_by_sweep(old_kt);
  else
    if (global.solver == ADR)
      kt = find_kt_from_scratch_by_ADR();
    else
      kt = find_kt_from_scratch_by_track();

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt_from_scratch_by_ADR
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::find_kt_from_scratch_by_ADR()
{
  // Determine reflection coefficients of walls.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  Complex R_left  = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_right = r_wall ? r_wall->get_R12() : -1.0;

  // Find min and max refractive index of structure.

  Complex min_eps_mu = materials[0]->eps_mu();
  Complex max_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;

    if (real(eps_mu) > real(min_eps_mu))
      max_eps_mu = eps_mu;
  }

  // Locate zeros starting from initial contour.

  const Real C = pow(2*pi/global.lambda, 2) / eps0 / mu0;
  
  Real max_kt = abs(sqrt(C*(max_eps_mu - min_eps_mu)))+1;

  Real Re = real(global.C_upperright);
  Real Im = imag(global.C_upperright);

  ExpandDirection dir;
  Complex lowerleft, upperright;
  if (real(min_eps_mu) > 0) // No metals.
  {
    lowerleft  = Complex(-0.1,      -0.1);
    upperright = Complex( Re*max_kt, Im*max_kt);
    dir = ur;
  }
  else
  {
    lowerleft  = Complex(-0.1,      -Im*max_kt);
    upperright = Complex( Re*max_kt, 0.1);
    dir = dr;
  }

  SlabDisp disp(materials, thicknesses, global.lambda, l_wall, r_wall);
  unsigned int zeros = global.N + 2 + materials.size();
  vector<Complex> kt = N_roots(disp, zeros, lowerleft, upperright,
                               1e-4, 1e-4, 4, dir);
  
  //cout << "Calls to slab dispersion relation : "<<disp.times_called()<<endl;
  
  // Eliminate false zeros.

  const Real eps_copies = 1e-6;
  
  for (unsigned int i=0; i<materials.size(); i++)
  {
    Complex kt_i = sqrt(C*(materials[i]->eps_mu() - min_eps_mu));
    
    remove_elems(&kt,     kt_i,      eps_copies);
    remove_elems(&kt,    -kt_i,      eps_copies);
    remove_elems(&kt, Complex(0.0),  eps_copies);
  }
  
  // In case of a homogeneous medium, add a TEM mode if appropriate.
  
  if (materials.size() == 1)
    if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
                                       && (abs(R_right + 1.0) < 1e-10) )
      || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
                                       && (abs(R_right - 1.0) < 1e-10) ) )
      kt.insert(kt.begin(), 0);

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt_from_scratch_by_track
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::find_kt_from_scratch_by_track()
{  
  // Set constants.
  
  const Real eps        = 1e-13;
  const Real eps_zero   = 1e-10;
  const Real eps_copies = 1e-6;
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;
  
  // Determine reflection coefficients of walls.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  Complex R_left  = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_right = r_wall ? r_wall->get_R12() : -1.0;

  bool branchcut = false;
  if ( (abs(R_left) < eps_zero) || (abs(R_right) < eps_zero) )
    branchcut = true; // Branchcuts exist only in open structures.

  bool TBC = false; // Transparent boundary conditions.
  if (   dynamic_cast<SlabWall_TBC*>(l_wall)
      || dynamic_cast<SlabWall_TBC*>(r_wall) )
    TBC = true;

  bool PC = false; // Photonic crystal boundary conditions.
  if (   dynamic_cast<SlabWall_PC*>(l_wall)
      || dynamic_cast<SlabWall_PC*>(r_wall) )
    PC = true;

  // Find min and max refractive index of lossless structure.

  Real min_eps_mu_lossless = real(materials[0]->eps_mu());
  Real max_eps_mu_lossless = real(materials[0]->eps_mu());
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Real eps_mu_lossless = real(materials[i]->eps_mu());
    
    if (eps_mu_lossless < min_eps_mu_lossless)
      min_eps_mu_lossless = eps_mu_lossless;

    if (eps_mu_lossless > min_eps_mu_lossless)
      max_eps_mu_lossless = eps_mu_lossless;
  }
  
  // Create dispersion relation for lossless structure.
  
  SlabDisp disp(materials, thicknesses, lambda, l_wall, r_wall);
  params = disp.get_params();

  vector<Complex> params_lossless;
  for (unsigned int i=0; i<params.size(); i++)
    params_lossless.push_back(real(params[i]));

  disp.set_params(params_lossless);

  // Make a rough estimate of separation between modes.

  Real d_kt = real(pi/get_width());

  // Find the modes of the structure, proceeding in a number
  // of possible steps:
  //
  //   I:   find modes along imag axis for lossless structure.
  //   II:  find modes along real axis for lossless structrue.
  //   III: track those to modes of true lossy structure.
  //   IV:  find modes in complex plane for true structure.


  
  //
  // I: Find propagating modes of lossless structure.
  //
  
  Wrap_imag_to_abs prop_wrap(disp);

  const Real C = pow(2*pi/global.lambda, 2) / eps0 / mu0;  
  Real prop_kt_end_lossless 
    = abs(sqrt(C*(max_eps_mu_lossless - min_eps_mu_lossless)));

  vector<Real> kt_prop_lossless;
  if (abs(prop_kt_end_lossless) > 0)
  {
    if (d_kt > prop_kt_end_lossless)
      d_kt = prop_kt_end_lossless;

    kt_prop_lossless = brent_all_minima
      (prop_wrap,0.0001,prop_kt_end_lossless,d_kt/global.precision,eps,1);
  }
  else
    prop_kt_end_lossless = 0.25; // To make rectangle for complex zero search.
  
  std::reverse(kt_prop_lossless.begin(), kt_prop_lossless.end());

  // Eliminate false zeros.

  remove_elems(&kt_prop_lossless, 0.0, eps_copies);
  for (unsigned int i=0; i<materials.size(); i++)
  {
    Real kt_i = sqrt(abs(C*(real(materials[i]->eps_mu())
                               - min_eps_mu_lossless)));

    remove_elems(&kt_prop_lossless,  kt_i, eps_copies);
    remove_elems(&kt_prop_lossless, -kt_i, eps_copies);
  }
  
  // Check if the minima found correspond really to zeros by using a mueller
  // Don't do this when there are branchcuts, as this can cause problems for
  // numerical stability.

  vector<Complex> kt_lossless;
  for (unsigned int i=0; i<kt_prop_lossless.size(); i++)
  { 
    if ( branchcut && (abs(disp(I*kt_prop_lossless[i])) > 1e-2) )
    {
      std::ostringstream s;
      s << "Warning: possibly insufficient precision around kt "
        << kt_prop_lossless[i] << "." << std::endl;
      s << "Removing this candidate.";
      py_print(s.str());
    }
    else
    {
      bool error = false;

      Real kt_new = imag(mueller(disp, I*kt_prop_lossless[i],
                                 I*kt_prop_lossless[i]+0.002,1e-14,
                                 0,100,&error));      
      if (!error)
        kt_lossless.push_back(I*kt_new);
    }
  }



  //
  // II: Find evanescent modes of lossless structure.
  //

  if (!(branchcut || TBC))
  {
    Wrap_real_to_abs evan_wrap(disp);
    
    Real kt_begin = .0001;
    int extra = 1;
    int modes_left = (kt_lossless.size() >= global.N + 2) ? 0
      : global.N - kt_lossless.size() + 2;
    
    while (modes_left)
    {
      if (extra > 1)
        py_print("Lost too many candidate modes. Expanding search region.");
      
      vector<Real> kt_evan_lossless = brent_N_minima
        (evan_wrap,kt_begin,modes_left+extra,d_kt/global.precision_rad,eps,1);
      
      // Check if the minima found correspond really to zeros by using 
      // a mueller solver.

      for (unsigned int i=0; i<kt_evan_lossless.size(); i++)
      {
        bool error = false;

        Real kt_new = branchcut
          ? kt_evan_lossless[i]
          : real(mueller(disp,kt_evan_lossless[i],
                         kt_evan_lossless[i]+0.002*I,1e-14,0,100,&error));

        if (!error)
          kt_lossless.push_back(kt_new);
      }

      modes_left = (kt_lossless.size() >= global.N + 2) ? 0
        : global.N - kt_lossless.size() + 2;

      kt_begin = real(kt_lossless[kt_lossless.size()-1])+.0001;
      extra++;
    }
  }

  
  
  //
  // III: Eliminate double and false zeros, find degenerate zeros
  //      snap them to axis and trace modes to those of the true structure.
  //

  remove_copies   (&kt_lossless, eps_copies);
  remove_opposites(&kt_lossless, eps_copies);
  remove_elems    (&kt_lossless, Complex(0.0), eps_copies);

  for (unsigned int i=0; i<materials.size(); i++)
  {
    Complex kt_i = sqrt(C*(real(materials[i]->eps_mu())-min_eps_mu_lossless));

    remove_elems(&kt_lossless,  kt_i, eps_copies);
    remove_elems(&kt_lossless, -kt_i, eps_copies);
  }

  vector<Complex> kt_lossless_single = kt_lossless;
  if (global.degenerate)
    kt_lossless = mueller_multiple(disp, kt_lossless_single, 1e-12, 0, 100);
  for (unsigned int i=0; i<kt_lossless.size(); i++)
  {
    if (abs(real(kt_lossless[i])) < abs(imag(kt_lossless[i])))
      kt_lossless[i] = I*imag(kt_lossless[i]);
    else
      kt_lossless[i] =   real(kt_lossless[i]);
  }

  bool degenerate = kt_lossless.size() > kt_lossless_single.size();
  
  vector<Complex> kt, forbidden;
  if (global.chunk_tracing && !degenerate)
    kt = traceroot_chunks
      (kt_lossless,disp,params_lossless,params,forbidden,global.sweep_steps);
  else
    kt = traceroot
      (kt_lossless,disp,params_lossless,params,forbidden,global.sweep_steps);



  //
  // IV: Find modes in complex plane
  //

  disp.set_params(params);

  if ( (branchcut || TBC) && (kt.size() < global.N) )
  {
    Real evan_kt_end = 2*(global.N-kt.size())*d_kt; // Not general!

    Complex lowerleft(0.001, 0.001);
    Complex upperright = global.C_upperright;

    const int sections = int(global.C_steps);

    vector<Complex> kt_complex = allroots(disp, lowerleft, upperright);

    disp.set_params(params);

    // Eliminate doubles and false zeros.

    remove_copies(&kt_complex, eps_copies);

    for (unsigned int i=0; i<materials.size(); i++)
    {
      Complex kt_i = sqrt(C*(materials[i]->eps_mu() - disp.get_min_eps_mu()));

      remove_elems(&kt_complex,  kt_i, eps_copies);
      remove_elems(&kt_complex, -kt_i, eps_copies);
    }

    std::ostringstream s;
    s << "Found " << kt_complex.size() << " complex zeros in region "
      << lowerleft << " " << upperright << ".";
    py_print(s.str());

    for (unsigned int i=0; i<kt_complex.size(); i++)
      kt.push_back(kt_complex[i]);
  }



  //
  // Finalise.
  //
  
  // In case of a homogeneous medium, add a TEM mode if appropriate.
  
  if (materials.size() == 1)
    if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
                                       && (abs(R_right + 1.0) < 1e-10) )
      || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
                                       && (abs(R_right - 1.0) < 1e-10) ) )
      kt.insert(kt.begin(), 0);

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt_by_sweep
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::find_kt_by_sweep(vector<Complex>& old_kt)
{ 
  // Set constants.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;
    
  // Trace modes from old configuration to new one.
  
  SlabDisp disp(materials, thicknesses, global.lambda, l_wall, r_wall);
  vector<Complex> params_new = disp.get_params();
    
  vector<Complex> forbidden;
  vector<Complex> kt =
    traceroot(old_kt, disp, params, params_new, forbidden, global.sweep_steps);

  params = params_new;

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::build_modeset
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::build_modeset(const vector<Complex>& kt)
{
  // Check if enough modes are found.

  if (kt.size() < global.N)
  {
    std::ostringstream s;
    s << "Error: didn't find enough modes ("
      << kt.size() << "/" << global.N << "). ";
    py_error(s.str());
    exit (-1);
  }
  
  // Clear old modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Find minimum eps mu.

  Complex min_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }
  
  const Real C = pow(2*pi/global.lambda, 2) / (eps0 * mu0);

  // Create new modeset.
  
  bool calc_fw = true;
  for (unsigned int i=0; i<kt.size(); i++)
  {
    Complex kz = sqrt(C*min_eps_mu - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    Polarisation pol = global.polarisation;
    if (global.polarisation == TE_TM)
      pol = ( i < int(kt.size()/2) ) ? TE : TM;
    
    calc_fw = !calc_fw;
    Slab_M_Mode *newmode = new Slab_M_Mode(pol, kz, kt[i], this, calc_fw);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat; 
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::get_params() const
{
  vector<Complex> params;

  for (unsigned int i=0; i<thicknesses.size(); i++)
  {
    params.push_back(thicknesses[i]);
    params.push_back(materials[i]->n());
    params.push_back(materials[i]->mur());
  }
  
  return params;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::set_params
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::set_params(const vector<Complex>& params)
{
  unsigned int params_index = 0;
  Complex current_x = 0.0;
    
  for (unsigned int i=0; i<thicknesses.size(); i++)
  {
    thicknesses[i] = params[params_index++];
    
    current_x += thicknesses[i];
    discontinuities[i] = current_x;
    
    materials[i]->set_n  (params[params_index++]);
    materials[i]->set_mur(params[params_index++]);
  }

  last_lambda = 0;
  slabmatrix_cache.deregister(this);

  // Determine core.

  unsigned int core_index = 0;
  for (unsigned int i=0; i<materials.size(); i++)
    if ( real(materials[i]->eps_mu()) > real(materials[core_index]->eps_mu()) )
      core_index = i;

  core = materials[core_index];
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlab::find_modes()
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

  // Only TE or TM modes needed.

  if ((global.polarisation == TE) || (global.polarisation == TM))
  {
    vector<Complex> kt(find_kt());
    build_modeset(kt);  
  }

  // Both TE and TM modes needed.

  if (global.polarisation == TE_TM)
  {
    // Cheat on global variables.

    const int n = int(global.N/2);
    if (2*n != global.N)
      py_print("Warning: changing N to even number.");
    global.N = n;

    const Complex old_beta = global_slab.beta;
    global_slab.beta = 0.0;

    // Find TE modes.

    global.polarisation = TE;
    vector<Complex> kt(find_kt());

    // Find TM modes

    global.polarisation = TM;
    vector<Complex> kt_TM(find_kt());

    // Restore global variables and build modeset.

    global.N = 2*n;
    global_slab.beta = old_beta;
    global.polarisation = TE_TM;

    kt.insert(kt.end(), kt_TM.begin(), kt_TM.end());
    build_modeset(kt);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::find_kt
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> UniformSlab::find_kt()
{
  vector<Complex> kt;

  // Determine reflection coefficients of walls.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  Complex R_left  = l_wall ? l_wall->get_R12() : -1;
  Complex R_right = r_wall ? r_wall->get_R12() : -1;

  // Check values.

  bool analytic = true;
  if ( ((abs(R_left -1.0) > 1e-10) && (abs(R_left +1.0) > 1e-10)) ||
       ((abs(R_right-1.0) > 1e-10) && (abs(R_right+1.0) > 1e-10)) )
    analytic = false;
  
  if (global.lambda == 0)
  {
    py_error("Error: wavelength not set.");
    return kt;
  }
  
  if (global.N == 0)
  {
    py_error("Error: number of modes not set.");
    return kt;
  }

  // If no analytic solution is available, follow same route as
  // non-uniform slab.

  if (!analytic)
  {
    Slab_M tmp = (*core)(get_width()/2.) + (*core)(get_width()/2.);

    if (l_wall)
      tmp.set_left_wall(*l_wall);

    if (r_wall)
      tmp.set_right_wall(*r_wall);
   
    vector<Complex> old_kt; // Sweep not implemented for this case.
    return tmp.find_kt(old_kt);
  }

  // Use analytical solution if available.

  // TEM mode.
  
  if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
                                     && (abs(R_right + 1.0) < 1e-10) )
    || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
                                     && (abs(R_right - 1.0) < 1e-10) ) )
    kt.push_back(0.0);

  // Other modes.
  
  Complex start  = (abs(R_left*R_right-1.0) < 1e-10) ? 0 : pi;
  unsigned int i = (abs(R_left*R_right-1.0) < 1e-10) ? 1 : 0;

  while (kt.size() != global.N)
  {
    kt.push_back((start+2*i*pi)/2./get_width());
    i++;
  }

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::build_modeset
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlab::build_modeset(vector<Complex>& kt)
{  
  // Clear old modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Create new modeset.

  const Complex k = 2*pi/global.lambda * core->n();
  
  for (unsigned int i=0; i<kt.size(); i++)
  { 
    Complex kz = sqrt(k*k - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    Polarisation pol = global.polarisation;
    if (global.polarisation == TE_TM)
      pol = ( i < int(kt.size()/2) ) ? TE : TM;

    UniformSlabMode *newmode = new UniformSlabMode(pol, kz, this);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat; 
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::get_params
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> UniformSlab::get_params() const
{
  vector<Complex> params;

  params.push_back(discontinuities[0]);
  params.push_back(core->n());
  params.push_back(core->mur());

  return params;
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::set_params
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlab::set_params(const vector<Complex>& params)
{
  discontinuities[0] = params[0];
  core->set_n(params[1]);
  core->set_mur(params[2]);

  last_lambda = 0;
  slabmatrix_cache.deregister(this);
}
