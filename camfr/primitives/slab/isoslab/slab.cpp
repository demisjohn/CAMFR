
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
#include <fstream>
#include "slab.h"
#include "slabwall.h"
#include "slabdisp.h"
#include "slabmode.h"
#include "slaboverlap.h"
#include "../slabmatrixcache.h"
#include "../../planar/planar.h"
#include "../../../math/calculus/calculus.h"
#include "../../../math/calculus/fourier/fourier.h"

using std::vector;

#include "../../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::Slab_M
//
/////////////////////////////////////////////////////////////////////////////

Slab_M::Slab_M(const Expression& expression, int M_series_) 
  : M_series(M_series_)
{ 
  Expression ex = expression.flatten();
  Complex current_x = 0.0;

  for (unsigned int i=0; i<ex.get_size(); i++)
  {
    Material* m = dynamic_cast<Material*>(ex.get_term(i)->get_mat());

    if (!m)
    {
      py_error("Error: expression contains non-material term.");
      return;
    }

    Complex thickness = ex.get_term(i)->get_d();

    if (i == 0)
      thickness += I*global_slab.lower_PML;

    if (i == ex.get_size()-1)
      thickness += I*global_slab.upper_PML;

    // Combine succesive terms containing the same material.      

    while ( (i+1 < ex.get_size())
         && (*m == *(dynamic_cast<Material*>(ex.get_term(i+1)->get_mat()))) )
    {
      thickness += ex.get_term(i+1)->get_d();
    
      if (i+1 == 0)
        thickness += I*global_slab.lower_PML;

      if (i+1 == ex.get_size()-1)
        thickness += I*global_slab.upper_PML;
      
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

  // Check if waveguide is uniform.

  uniform = true;
  for (unsigned int i=0; i<materials.size(); i++)    
    if (    (abs(materials[i]->n()   - core->n())   > 1e-12)
         || (abs(materials[i]->mur() - core->mur()) > 1e-12) )
      uniform = false;

  //SlabWall* l_wall_disp = lowerwall ? lowerwall : global_slab.lowerwall;
  //SlabWall* u_wall_disp = upperwall ? upperwall : global_slab.upperwall;
  //disp_f = new SlabDisp(materials,thicknesses, global.lambda, 
  //                      l_wall_disp,u_wall_disp);

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
// Slab_M::material_at
//
/////////////////////////////////////////////////////////////////////////////

Material* Slab_M::material_at(const Coord& coord) const
{
  unsigned int i = index_lookup(coord.c1, coord.c1_limit, discontinuities);
  
  if (i == materials.size())
  {
    py_error("Warning: coordinate out of range in eps_at. Restricting it.");
    i = materials.size()-1;
  }

  return materials[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::eps_avg
//
/////////////////////////////////////////////////////////////////////////////

Complex Slab_M::eps_avg() const
{
  Complex avg = 0.0;

  for (unsigned int i=0; i<materials.size(); i++)
    avg += materials[i]->eps() * thicknesses[i];
  
  return avg / get_width();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::is_mirror_image_of
//
/////////////////////////////////////////////////////////////////////////////

bool Slab_M::is_mirror_image_of(const SlabImpl* m) const
{
  const Slab_M* medium_II = dynamic_cast<const Slab_M*>(m);
  if (!medium_II)
    return false;

  const unsigned int K = thicknesses.size();
  if (K != medium_II->thicknesses.size())
    return false;

  const Real eps = 1e-8;

  for (int i=0; i<K; i++)
    if (    (abs(thicknesses[i] - medium_II->thicknesses[K-1-i]) > eps)
         || (*materials[i] != *(medium_II->materials[K-1-i])) )
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes()
{ 
  // Check values.

  if (real(global.lambda) == 0)
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

    const Complex old_beta = global.slab_ky;
    global.slab_ky = 0.0;

    vector<Complex> old_params = params;

    // Find TE modes.

    global.polarisation = TE;

    vector<Complex> old_kt_TE;
    if (    modeset.size() 
         && (modeset.front()->pol == TE) 
         && (modeset.back() ->pol == TM) )
      for (unsigned int i=0; i<n; i++)
        old_kt_TE.push_back(dynamic_cast<Slab_M_Mode*>(modeset[i])->get_kt());

    vector<Complex> kt(find_kt(old_kt_TE));
    kt.erase(kt.begin()+n, kt.end());

    // Find TM modes

    params = old_params;

    last_lambda = 0.0; // Force a recalc.

    global.polarisation = TM;

    vector<Complex> old_kt_TM;
    if (    modeset.size() 
         && (modeset.front()->pol == TE) 
         && (modeset.back() ->pol == TM) )
      for (unsigned int i=n; i<2*n; i++)
        old_kt_TM.push_back(dynamic_cast<Slab_M_Mode*>(modeset[i])->get_kt());

    vector<Complex> kt_TM(find_kt(old_kt_TM));

    // Restore global variables and build modeset.

    global.N = 2*n;
    global.slab_ky = old_beta;
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
  // from scratch. Note that the series solvers are currently only 
  // fully implemented for TM.

  if (     (global.solver == series) 
        || (global.solver == ASR            && global.polarisation == TM) 
        || (global.solver == stretched_ASR  && global.polarisation == TM) )
      return find_kt_from_estimates();

  if (     (global.solver == ASR            && global.polarisation == TE)
        || (global.solver == stretched_ASR  && global.polarisation == TE) )
      return find_kt_from_scratch_by_track();
     
  if (global.sweep_from_previous && (modeset.size() >= global.N))
      return find_kt_by_sweep(old_kt);
  
  if (global.solver == ADR)
    return find_kt_from_scratch_by_ADR();

  return find_kt_from_scratch_by_track();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt_from_scratch_by_ADR
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Slab_M::find_kt_from_scratch_by_ADR()
{
  // Determine reflection coefficients of walls.

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;

  // Find min and max refractive index of structure.

  Complex min_eps_mu = materials[0]->eps_mu();
  Complex max_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;

    if (real(eps_mu) > real(max_eps_mu))
      max_eps_mu = eps_mu;
  }

  // Locate zeros starting from initial contour.

  const Complex C = pow(2*pi/global.lambda, 2) / eps0 / mu0;
  
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

  SlabDisp disp(materials,thicknesses,global.lambda,l_wall,u_wall);
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
  
  if (uniform)
    if ( ( (global.polarisation == TM) && (abs(R_lower + 1.0) < 1e-10)
                                       && (abs(R_upper + 1.0) < 1e-10) )
      || ( (global.polarisation == TE) && (abs(R_lower - 1.0) < 1e-10)
                                       && (abs(R_upper - 1.0) < 1e-10) ) )
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
  
  const Complex lambda  = global.lambda;
  const Complex k0      = 2*pi/lambda;
  
  // Determine reflection coefficients of walls.

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;

  bool branchcut = false;
  if ( (abs(R_lower) < eps_zero) || (abs(R_upper) < eps_zero) )
    branchcut = true; // Branchcuts exist only in open structures.

  bool TBC = false; // Transparent boundary conditions.
  if (   dynamic_cast<SlabWall_TBC*>(l_wall)
      || dynamic_cast<SlabWall_TBC*>(u_wall) )
    TBC = true;

  bool PC = false; // Photonic crystal boundary conditions.
  if (   dynamic_cast<SlabWall_PC*>(l_wall)
      || dynamic_cast<SlabWall_PC*>(u_wall) )
    PC = true;

  // Find min and max refractive index of lossless structure.

  Real min_eps_mu_lossless = real(materials[0]->eps_mu());
  Real max_eps_mu_lossless = real(materials[0]->eps_mu());
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Real eps_mu_lossless = real(materials[i]->eps_mu());
    
    if (eps_mu_lossless < min_eps_mu_lossless)
      min_eps_mu_lossless = eps_mu_lossless;

    if (eps_mu_lossless > max_eps_mu_lossless)
      max_eps_mu_lossless = eps_mu_lossless;
  }

  bool metallic = (min_eps_mu_lossless < 0.0);
  
  // Create dispersion relation for lossless structure.

  SlabDisp disp(materials,thicknesses,lambda,l_wall,u_wall);
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

  const Complex C = pow(2*pi/global.lambda, 2) / eps0 / mu0;

  Complex max_eps_eff = max_eps_mu_lossless;

  if (    (global.polarisation == TM) // Surface plasmon.
       && (max_eps_mu_lossless*min_eps_mu_lossless < 0.0) )
    max_eps_eff = 1.0/(1.0/max_eps_mu_lossless + 1.0/min_eps_mu_lossless);

  Real prop_kt_end_lossless = global_slab.estimate_cutoff 
                         *abs(sqrt(C*(max_eps_eff - min_eps_mu_lossless)));
  
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

  //for (int i=0; i<kt_prop_lossless.size(); i++)
  //  std::cout << i << " " << kt_prop_lossless[i] << std::endl;

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

  vector<Complex> kt_lossless, kt_complex;
  for (unsigned int i=0; i<kt_prop_lossless.size(); i++)
  {
    const Real f = abs(disp(I*kt_prop_lossless[i]));
    
    if ( (f > 1.0) || (branchcut && (f > 1e-2)) )
    {
      std::ostringstream s;
      s << "Warning: possibly insufficient precision around kt "
        << kt_prop_lossless[i] << "." << std::endl;
      s << "Removing this candidate with f(x)=" << f << ".";
      py_print(s.str());
    }
    else
    {
      bool error = false;

      Complex kt_new = mueller(disp, I*kt_prop_lossless[i],
                               I*kt_prop_lossless[i]+0.002,1e-11,0,100,&error);

      if (!error && metallic && abs(real(kt_new)) > 0.001)
      {
        py_print("Found complex mode pair.");
        kt_complex.push_back(kt_new);
      }
      else if (!error)
        kt_lossless.push_back(I*imag(kt_new));
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

        Complex kt_new = branchcut
          ? kt_evan_lossless[i]
          : mueller(disp,kt_evan_lossless[i],kt_evan_lossless[i]+0.002*I,
                    1e-11,0,100,&error);

        if (!error && metallic && (abs(imag(kt_new)) > 0.001))
        {
          py_print("Found complex mode pair.");
          kt_complex.push_back(kt_new);
        }
        else if (!error)
          kt_lossless.push_back(real(kt_new));
      }

      modes_left = (kt_lossless.size() + 2*kt_complex.size() >= global.N + 2) 
        ? 0
        : global.N - kt_lossless.size() -2*kt_complex.size() + 2;

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
    kt_lossless = mueller_multiple(disp, kt_lossless_single, 1e-11, 0, 100);
  for (unsigned int i=0; i<kt_lossless.size(); i++)
  {
    if (abs(real(kt_lossless[i])) < abs(imag(kt_lossless[i])))
      kt_lossless[i] = I*imag(kt_lossless[i]);
    else
      kt_lossless[i] =   real(kt_lossless[i]);
  }

  remove_elems(&kt_lossless, Complex(0.0), eps_copies);
  
  bool degenerate = kt_lossless.size() > kt_lossless_single.size();

  for (unsigned int i=0; i<kt_complex.size(); i++)
  {  
    kt_lossless.push_back(kt_complex[i]);
    kt_lossless.push_back(-real(kt_complex[i])+I*imag(kt_complex[i]));
  }

  vector<Complex> kt, forbidden;
  forbidden.push_back(Complex(0.0));
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
      Complex kt_i = sqrt(C*(materials[i]->eps_mu() - disp.get_kt_eps_mu()));

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
  
  if (uniform)
    if ( ( (global.polarisation == TM) && (abs(R_lower + 1.0) < 1e-10)
                                       && (abs(R_upper + 1.0) < 1e-10) )
      || ( (global.polarisation == TE) && (abs(R_lower - 1.0) < 1e-10)
                                       && (abs(R_upper - 1.0) < 1e-10) ) )
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

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;
    
  // Trace modes from old configuration to new one.
  
  SlabDisp disp(materials, thicknesses, global.lambda, l_wall, u_wall);
  vector<Complex> params_new = disp.get_params();
    
  vector<Complex> forbidden;
  vector<Complex> kt =
    traceroot(old_kt, disp, params, params_new, forbidden, global.sweep_steps);

  params = params_new;

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// fill_E_matrix
//
//  Create matrix for a TM RCWA eigenvalue problem.
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::fill_E_matrix(cMatrix* E, int M, int n, 
                           const Complex& lambda_over_d, 
                           const Complex& offset) const
{
  cVector f_eps    (4*M+1,fortranArray);
  cVector f_inv_eps(4*M+1,fortranArray);

  // Construct data structures for fourier analysis.

  vector<Complex> disc(discontinuities);
  disc.insert(disc.begin(), 0.0);

  vector<Complex> eps;
  vector<Complex> h_eps;
  vector<Complex> inv_eps;

  for (int i=0; i<disc.size()-1; i++)
  {
    Complex eps_i = eps_at(Coord(disc[i], 0, 0, Plus))/eps0;
 
    eps.push_back(eps_i);
    inv_eps.push_back(1.0/eps_i);
    h_eps.push_back(1.0);
  }

  // Construct fourier matrices.
 
  bool extend = true;
  bool stretch = (global.solver == stretched_ASR);
  Real eta = global_slab.eta_ASR;

  if ( (global.solver == ASR) || (global.solver == stretched_ASR) )
  {
   
   f_eps     = fourier_ASR(    eps, disc, 2*M, 0, stretch, extend, eta);
   f_inv_eps = fourier_ASR(inv_eps, disc, 2*M, 0, stretch, extend, eta);
  }
  else
  {
    f_eps     = fourier(    eps, disc, 2*M, 0, extend);
    f_inv_eps = fourier(inv_eps, disc, 2*M, 0, extend);
  }

  cMatrix Eps(n,n,fortranArray);
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      Eps(i,j) = f_eps(i-j + 2*M+1);

  cMatrix A(n,n,fortranArray);
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      A(i,j) = f_inv_eps(i-j + 2*M+1);

  // Note: the following inversions can be sped up because these
  // are Toeplitz matrices.

  cMatrix inv_Eps(n,n,fortranArray); inv_Eps.reference(invert(Eps));
  cMatrix inv_A  (n,n,fortranArray); inv_A  .reference(invert(A));

  for (int i=1; i<=n; i++)
  {
    for (int j=1; j<=n; j++)
      (*E)(i,j) = (Real(i-1-M)*lambda_over_d+offset) * inv_Eps(i,j) *
                  (Real(j-1-M)*lambda_over_d+offset);
  }

  // Subtract the unity matrix or the ASR version of the unity matrix.

  if ((global.solver != ASR) && (global.solver != stretched_ASR))
  {
    for (int i=1; i<=n; i++)
      (*E)(i,i) -= 1.0;
  }
  else
  {
    cVector f_h_eps(4*M+1,fortranArray);
    f_h_eps = fourier_ASR(h_eps, disc, 2*M, 0, extend, stretch, eta);

    for (int i=1; i<=n; i++)
      for (int j=1; j<=n; j++)
        (*E)(i,j) -= f_h_eps(i-j + 2*M+1);
  }
  
  *E = multiply(inv_A, *E);
}



/////////////////////////////////////////////////////////////////////////////
//
// fill_E_matrix_biaxial
//
//  Create matrix for TM RCWA problem eigenvalue with biaxial media.
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::fill_E_matrix_biaxial(cMatrix* E, int M, int n,
                                   const Complex& lambda_over_d,
                                   const Complex& offset) const
{
  // Calculate fourier transform of eps(x) profile.

  vector<Complex> disc_x(get_discontinuities());
  disc_x.insert(disc_x.begin(), 0.0);

  vector<Complex> eps, inv_eps, mu;
  for (int k=0; k<disc_x.size()-1; k++)
  {
    Material* mat = material_at(Coord(disc_x[k],0,0,Plus));

    eps.    push_back(    mat->epsr(3));
    inv_eps.push_back(1.0/mat->epsr(1));
    mu.     push_back(    mat->mur (2));
  }

  // TODO: check if it is an advantage to join this in a single function call.

  cVector f_eps    (4*M+1,fortranArray);
  cVector f_inv_eps(4*M+1,fortranArray);
  cVector f_mu     (4*M+1,fortranArray);

  f_eps     = fourier(eps,     disc_x, 2*M, NULL, true);
  f_inv_eps = fourier(inv_eps, disc_x, 2*M, NULL, true);
  f_mu      = fourier(mu,      disc_x, 2*M, NULL, true);

  cMatrix Eps(n,n,fortranArray);
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      Eps(i,j) = f_eps(i-j + 2*M+1);

  cMatrix A(n,n,fortranArray);
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      A(i,j) = f_inv_eps(i-j + 2*M+1);

  // Note: the following inversions can be sped up because these
  // are Toeplitz matrices.

  cMatrix inv_Eps(n,n,fortranArray); inv_Eps.reference(invert(Eps));
  cMatrix inv_A  (n,n,fortranArray); inv_A  .reference(invert(A));

  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      (*E)(i,j) = (Real(i-1-M)*lambda_over_d+offset) * inv_Eps(i,j) *
        (Real(j-1-M)*lambda_over_d+offset) - f_mu(i-j + 2*M+1);
  
  *E = multiply(inv_A, *E);
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::estimate_kz2_from_uniform_modes
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> Slab_M::estimate_kz2_from_uniform_modes()
{  
  // Set constants.

  const Complex omega = 2*pi/global.lambda * c;
  const Complex k0 = 2*pi/global.lambda;

  // Create eigenvalue problem.

  if (global.polarisation == TE) // TODO: reformulate as fourier expansion.
  {
    const int n = M_series ? M_series : int(global.N * global.mode_surplus);
  
    cMatrix E(n,n,fortranArray);

    // Calculate overlap matrices.

    int old_N = global.N;
    global.N = n;

    Material ref_mat(1.0);
    UniformSlab ref(get_width(),ref_mat);
    ref.find_modes();
   
    global.N = old_N;

    overlap_reference_modes(&E, NULL, ref, *this);
    
    for (int i=1; i<=n; i++)
    {
      SlabMode* TE_i = dynamic_cast<SlabMode*>(ref.get_mode(i));

      Complex kz_i = TE_i->get_kz();
      if (abs(kz_i) < 1e-6)
        py_print("WARNING: reference mode close to cutoff!");
      
      const Complex factor = omega * kz_i;
      for (int j=1; j<=n; j++)
        E(i,j) *= factor;

      E(i,i) -= pow(TE_i->kx_at(Coord(0,0,0)), 2);
    }

    // Solve eigenvalue problem.

    cVector e(n,fortranArray);
    if (global.stability == normal)
      e = eigenvalues(E);
    else
      e = eigenvalues_x(E);

    vector<Complex> kz2;
    for (int i=1; i<=n; i++)
      kz2.push_back(e(i));

    return kz2;

  }
  else // TM case.
  {
    // Determine reflection coefficients of walls.

    SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
    SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

    Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
    Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;

    Complex D = get_width();

    // Recalculate for a stretched structure.
    // TODO_PDB: this is probably where the udisc calculation should go?

    if (global.solver == stretched_ASR)
    {
      vector<Complex> disc = discontinuities;
      disc.insert(disc.begin(), 0.0);

      int S = disc.size();
      Complex begin = (S-1)*real(disc[1]  -disc[0])   + imag(disc[S-1])*I;
      Complex end   = (S-1)*real(disc[S-1]-disc[S-2]) + imag(disc[S-1])*I;

      D = (real(begin) >= real(end)) ? begin : end;
    }

    const Complex lambda_over_d = global.lambda / D / 2.;


    Complex offset = (abs(R_lower*R_upper-1.0) < 1e-10) ? 0 : lambda_over_d/2.;

    // Calculate eps matrices.

    int M = M_series ? M_series : int(global.N * global.mode_surplus);

    if ( (abs(R_lower+1.0) < 1e-6) && (abs(R_upper+1.0) < 1e-6) ) // EE
      M--; // Note: don't do this for HH since it has a dummy solution.

    int n = 2*M + 1;

    if (abs(offset) > 1e-6) // EH or HE has no zeroth order term.
      n--;
 
    cMatrix E(n,n,fortranArray);
    
    //if (global.eigen_calc != arnoldi) // Hijacked switch.
    //{	std::cout<<"Fill_E_matrix_biaxial"<<std::endl;
    //  fill_E_matrix_biaxial(&E, M, n, lambda_over_d, offset);
    //}
    // else
    // {std::cout<<"Fill_E_matrix"<<std::endl;

    fill_E_matrix(&E, M, n, lambda_over_d, offset);

    //  }

    // Reduce eigenvalue problem to retain only solutions with desired
    // symmetry.
    //
    // For same walls on both sides, the 2M+1 eigenvectors come in 2 sets:
    //
    //   [...,b,a,0,-a,-b,...] (HH)
    //   [...,b,a,c, a, b,...] (EE)
    //
    // For different walls and 2M eigenvectors we have:
    //
    //   [...,b,a,-a,-b,...] (HE) 
    //   [...,b,a, a, b,...] (EH)

    // TODO: try doing this reduction before matrix inversions 
    // and multiplications.
    
    int M_ = int((n+1)/2);
    
    cMatrix E_(M_,M_,fortranArray);

    if (abs(offset) < 1e-6) // EE or HH
    {
      for (int i=1; i<=M_; i++)
      {
        for (int j=1; j<=M_-1; j++)
          E_(i,j) = E(i,j) - R_lower*E(i,E.columns()-j+1);

        E_(i,M_) = (abs(R_lower-1.0) < 1e-6) ? 0.0 : E(i,M_);
      }
    }
    else // EH or HE
    {
      for (int i=1; i<=M_; i++)
        for (int j=1; j<=M_; j++)
          E_(i,j) = E(i,j) - R_lower*E(i,E.columns()-j+1);
    }

    // Solve eigenvalue problem.

    cVector e(E_.rows(),fortranArray);
    if (global.stability == normal)
      e = eigenvalues(E_);
    else
      e = eigenvalues_x(E_);

    vector<Complex> kz2;
    for (int i=1; i<=e.rows(); i++)
      if (abs(e(i)) > 1e-6) // Zero is dummy solution for HH case.
        kz2.push_back(-k0*k0*e(i));

    return kz2;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Auxiliary function objects
//
/////////////////////////////////////////////////////////////////////////////

struct kz_sorter
{
    bool operator()(const Complex& kz_a, const Complex& kz_b)
    {return ( real(kz_a*kz_a) > real(kz_b*kz_b) );}
};

struct kz2_sorter
{
    bool operator()(const Complex& kz2_a, const Complex& kz2_b)
      {return ( real(kz2_a) > real(kz2_b) );}
};

struct kt_sorter
{   bool operator()(const Complex& kt_a, const Complex& kt_b)
      {return ( real(kt_a*kt_a) < real(kt_b*kt_b) );}
};

struct kt_to_neff : ComplexFunction
{
  Complex C;

  kt_to_neff(const Complex& C_) : C(C_) {}

  Complex operator()(const Complex& kt)
    {return sqrt(C - kt*kt)/2./pi*global.lambda;}
};



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_kt_from_estimates
//
/////////////////////////////////////////////////////////////////////////////

std::vector<Complex> Slab_M::find_kt_from_estimates()
{
  vector<Complex> kz2;

  if (user_kz2_estimates.size() != 0)
    kz2 = user_kz2_estimates;
  else
    kz2 = estimate_kz2_from_uniform_modes();

  //else if (global.eigen_calc != arnoldi) // Hijacked switch.
  //  kz2 = estimate_kz2_from_uniform_modes();
  //else
  //  kz2 = estimate_kz2_from_RCWA();

  //for (unsigned int i=0; i<kz2.size(); i++)
  //  std::cout << "raw " << i << " " 
  //            << sqrt(kz2[i])/2./pi*global.lambda 
  //            << kz2[i] << std::endl;

  // Find min and max eps mu.

  Complex min_eps_mu = materials[0]->eps_mu();
  Complex max_eps_mu = materials[0]->eps_mu();

  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;

    if (real(eps_mu) > real(max_eps_mu))
      max_eps_mu = eps_mu;
  }

  //bool metallic = (real(min_eps_mu) < 0);
  bool metallic = global_slab.low_index_core;

  Complex kt_eps_mu = metallic ? max_eps_mu : min_eps_mu;
 
  const Complex C0 = pow(2*pi/global.lambda, 2) / eps0 / mu0;

  Complex max_eps_eff;
  max_eps_eff = max_eps_mu;

  if (    (global.polarisation == TM) // Surface plasmon.
       && (real(max_eps_mu)*real(min_eps_mu) < 0.0) )
    max_eps_eff = 1.0/(1.0/max_eps_mu + 1.0/min_eps_mu);

  Real max_kz = real(2.0*pi/global.lambda*sqrt(max_eps_eff/eps0/mu0)) 
                * global_slab.estimate_cutoff; // Safety factor.
  
 // Refine estimates.

  vector<Complex> kz2_coarse;

  for (unsigned int i=0; i<kz2.size(); i++)
  if ( global.keep_all_1D_estimates || (real(sqrt(kz2[i])) < max_kz) )
      kz2_coarse.push_back(kz2[i]);

  std::sort(kz2_coarse.begin(), kz2_coarse.end(), kz2_sorter());

  if (kz2_coarse.size() > global.N+5)
    kz2_coarse.erase(kz2_coarse.begin()+global.N+5, kz2_coarse.end());

  vector<Complex> kt_coarse;
  for (unsigned int i=0; i<kz2_coarse.size(); i++)
  {
    Complex kt = sqrt(C0*kt_eps_mu - kz2_coarse[i]);

    if (imag(kt) < 0)
      kt = -kt;

    if (abs(imag(kt)) < 1e-10)
      if (real(kt) > 0)
        kt = -kt;
    kt_coarse.push_back(kt);
  }  

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  kt_to_neff transform(C0*kt_eps_mu);
  SlabDisp disp(materials,thicknesses,global.lambda,l_wall,u_wall);

  vector<Complex> kt = mueller(disp, kt_coarse, 1e-11, 100, &transform, 0);

  // Eliminate false zeros.

  for (unsigned int i=0; i<materials.size(); i++)
  {
    Complex kt_i = sqrt(C0*(real(materials[i]->eps_mu())-kt_eps_mu));

    remove_elems(&kt,  kt_i, 1e-6);
    remove_elems(&kt, -kt_i, 1e-6);
  }

  return kt;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::build_modeset
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::build_modeset(vector<Complex>& kt)
{
  // Set polarisation.

  Polarisation old_polarisation = global.polarisation;

  // Clear old modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Find min and max eps mu.

  Complex min_eps_mu = materials[0]->eps_mu();
  Complex max_eps_mu = materials[0]->eps_mu();

  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;

    if (real(eps_mu) > real(max_eps_mu))
      max_eps_mu = eps_mu;
  }

  //bool metallic = (real(min_eps_mu) < 0);
  bool metallic = global_slab.low_index_core;

  Complex kt_eps_mu = metallic ? max_eps_mu : min_eps_mu;
 
  const Complex C = pow(2*pi/global.lambda, 2) / (eps0 * mu0);

  // Create Planars.

  vector<Planar> planars;
  for (unsigned int i=0; i<materials.size(); i++)
    planars.push_back(Planar(*materials[i]));
  
  // Determine walls.

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;

  MonoScatterer* l_boundary;
  if (abs(R_lower - (-1.0)) < 1e-9)
    l_boundary = new E_Wall_Mono(planars[0]);
  else if (abs(R_lower - 1.0) < 1e-9)
    l_boundary = new H_Wall_Mono(planars[0]);
  else
  {
    py_error("Unsupported wall type in Slab_M::build_modeset.");
    exit (-1);
  }

  MonoScatterer* u_boundary;
  if (abs(R_upper - (-1.0)) < 1e-9)
    u_boundary = new E_Wall_Mono(planars.back());
  else if (abs(R_upper - 1.0) < 1e-9)
    u_boundary = new H_Wall_Mono(planars.back());
  else
  {
    py_error("Unsupported wall type in Slab_M::build_modeset.");
    exit (-1);
  }

  // Find cores and create Planar stacks for layers to the left and right
  // of each core.

  vector<int>   core_index;
  vector<Stack> lower_stacks;
  vector<Stack> upper_stacks;
  
  for (unsigned int i=0; i<materials.size(); i++)
  {
    // Is this layer a core?
    
    bool is_core = false;

    if (materials.size() == 1)
      is_core = true;
    
    if (i == 0)
    {
      if (real(materials[i]->eps_mu()) > real(materials[i+1]->eps_mu()))
        is_core = true;
    }
    else if (i == materials.size()-1)    
    {
      if (real(materials[i]->eps_mu()) > real(materials[i-1]->eps_mu()))
        is_core = true;
    }
    else
    {   
      if (    (real(materials[i]->eps_mu()) > real(materials[i-1]->eps_mu()))
           && (real(materials[i]->eps_mu()) > real(materials[i+1]->eps_mu())) )
        is_core = true;
      if (    (real(materials[i]->eps_mu()) < 0.0) ) // Metal core.
         is_core = true;
    }


    // If so, create lower and upper stacks.

    if (is_core)
    {
      core_index.push_back(i);
      
      // Lower stack.

      Expression lower_ex;        
      lower_ex += Term(planars[i](thicknesses[i]/2.));
      for (int j=i-1; j>=0; j--)
        lower_ex += Term(planars[j](thicknesses[j]));
      lower_ex += Term(*l_boundary);
      lower_stacks.push_back(Stack(lower_ex));

      // Upper stack.

      Expression upper_ex;
      upper_ex += Term(planars[i](thicknesses[i]/2.)); 
      for (int j=i+1; j<materials.size(); j++)
        upper_ex += Term(planars[j](thicknesses[j]));
      upper_ex += Term(*u_boundary);
      upper_stacks.push_back(Stack(upper_ex));
      
      //std::cout << "lower" << lower_ex << std::endl;
      //std::cout << "upper" << upper_ex << std::endl;
    }
  }

  // Find candidate cores to calculate modes from.
  
  if (global.polarisation != TE_TM)
    std::sort(kt.begin(), kt.end(), kt_sorter());

  vector<vector<int> > best_cores; 
  for (unsigned int i=0; i<kt.size(); i++)
  {
          
    vector<Real> error_i;
    Real best_error;
    for (unsigned int j=0; j<core_index.size(); j++)
    {
      Complex kz = sqrt(C*kt_eps_mu - kt[i]*kt[i]);
    
      if (real(kz) < 0) 
        kz = -kz;

      if (abs(real(kz)) < abs(imag(kz))) // Prefer lossy modes.
        if (imag(kz) > 0)
          kz = -kz;

      if (old_polarisation == TE_TM)
        global.polarisation = ( i < int(global.N/2) ) ? TE : TM;
      Planar::set_kt(kz);

      lower_stacks[j].calcRT();
      upper_stacks[j].calcRT();

      Real error = abs(lower_stacks[j].R12(0,0)*upper_stacks[j].R12(0,0)-1.0);

      error_i.push_back(error);

      if (j==0)
        best_error = error;

      if (error < best_error)
        best_error = error;

      //std::cout << kz/2./pi*global.lambda 
      //          << " " << i << " " << j << " " << error 
      //          << std::endl << std::flush;
    }

    vector<int> best_cores_i;
    for (unsigned int j=0; j<core_index.size(); j++)
      if (abs(error_i[j] - best_error) < 1e-6)
        best_cores_i.push_back(j);

    //std::cout << "bestcores_i";
    //for (int k=0; k<best_cores_i.size(); k++)
    //  std::cout << " " << best_cores_i[k];
    //std::cout << std::endl << std::flush;

    best_cores.push_back(best_cores_i);
  }

  // Create new modeset.

  bool calc_fw = true; // To be removed, no longer used.
  int offset = 0;

  for (unsigned int i=0; i<kt.size(); i++)
  {
    // Pick best core from candidates.

    int core;

    // Heuristic for identical cores.

    if (global.field_calc_heuristic == identical)
    {
      if (best_cores[i].size() == 1)
      {
        core = best_cores[i][0];
        offset = 0;
      }
      else
      {
        if (offset >= best_cores[i].size())
          offset = 0;
      
        core = best_cores[i][offset++];
      }
    }

    // Heuristic for symmetric structures.

    else if (global.field_calc_heuristic == symmetric)
    {
      if (2*int(i/2) == i)
        core = best_cores[i].front();
      else
        core = best_cores[i].back();
    }
    
    //Complex neff = sqrt(C*kt_eps_mu - kt[i]*kt[i])/2./pi*global.lambda;
    //std::cout << "Mode " << i << " " << neff << " : ";
    //for (int j=0; j<best_cores[i].size(); j++)
    //  std::cout << best_cores[i][j] << " ";
    //std::cout << " : " << core << std::endl << std::flush;

    // Calculate kz.

    Complex kz = sqrt(C*kt_eps_mu - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < abs(imag(kz))) // Prefer lossy modes.
      if (imag(kz) > 0)
        kz = -kz;

    // Create mode.

    Polarisation pol = old_polarisation;
    if (old_polarisation == TE_TM)
      pol = ( i < int(global.N/2) ) ? TE : TM;
    global.polarisation = pol;

    calc_fw = !calc_fw;

    Slab_M_Mode *newmode;

    if (global.eigen_calc == arnoldi) // Old method.
      newmode = new Slab_M_Mode(pol, kz, kt[i], this, calc_fw);
    else
      newmode = new Slab_M_Mode(pol, kz, kt[i], this,
                                lower_stacks[core], upper_stacks[core]);

    //std::cout << "Mode " << i << kt[i] << kz/2./pi*global.lambda
    //          << " " << pol << " " << calc_fw << std::endl;
    //std::cout<<"lower_stacks[core]"<<" "<<core<<std::endl;
    //std::cout<<"upper_stacks[core]"<<" "<<core<<std::endl;
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }

  sort_modes();
  truncate_N_modes();

  // Gracefully handle case when not enough modes were found.

  if (modeset.size() < global.N)
  {
    std::ostringstream s;
    s << "Error: didn't find enough modes ("
      << kt.size() << "/" << global.N << "). ";
    py_error(s.str());
    
    while (modeset.size() != global.N)
    {
      Slab_M_Mode *dummy = new Slab_M_Mode(unknown, 0.0, 0.0, this);
      modeset.push_back(dummy);
    }
  }

  // Clean up.

  delete l_boundary;
  delete u_boundary;

  global.polarisation = old_polarisation;

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat_n = global.gain_mat->n();
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

  if (real(global.lambda) == 0)
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

    const Complex old_beta = global.slab_ky;
    global.slab_ky = 0.0;

    // Find TE modes.

    global.polarisation = TE;
    vector<Complex> kt(find_kt());

    // Find TM modes

    global.polarisation = TM;
    vector<Complex> kt_TM(find_kt());

    // Restore global variables and build modeset.

    global.N = 2*n;
    global.slab_ky = old_beta;
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

  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;

  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;

  // Check values.

  bool analytic = true;
  if ( ((abs(R_lower-1.0) > 1e-10) && (abs(R_lower+1.0) > 1e-10)) ||
       ((abs(R_upper-1.0) > 1e-10) && (abs(R_upper+1.0) > 1e-10)) )
    analytic = false;
  
  if (real(global.lambda) == 0)
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
      tmp.set_lower_wall(*l_wall);

    if (u_wall)
      tmp.set_upper_wall(*u_wall);
   
    vector<Complex> old_kt; // Sweep not implemented for this case.
    return tmp.find_kt(old_kt);
  }

  // Use analytical solution if available.

  // TEM mode.
  
  if ( ( (global.polarisation == TM) && (abs(R_lower + 1.0) < 1e-10)
                                     && (abs(R_upper + 1.0) < 1e-10) )
    || ( (global.polarisation == TE) && (abs(R_lower - 1.0) < 1e-10)
                                     && (abs(R_upper - 1.0) < 1e-10) ) )
    kt.push_back(0.0);

  // Other modes.
  
  Complex start  = (abs(R_lower*R_upper-1.0) < 1e-10) ? 0 : pi;
  unsigned int i = (abs(R_lower*R_upper-1.0) < 1e-10) ? 1 : 0;

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
    last_gain_mat_n = global.gain_mat->n(); 
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
