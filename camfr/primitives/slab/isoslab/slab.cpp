
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

#include "slab.h"
#include "slabwall.h"
#include "slabdisp.h"
#include "slabmode.h"
#include "slaboverlap.h"
#include "../slabmatrixcache.h"
#include "../../../math/calculus/calculus.h"
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
      cerr << "Error: expression contains non-material term." << endl;
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
    if ( real(materials[i]->n()) > real(materials[core_index]->n()) )
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
// Slab_M::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes()
{
  // Check values.

  if (global.lambda == 0)
  {
    cout << "Error: wavelength not set." << endl;
    return;
  }
  
  if (global.N == 0)
  {
    cout << "Error: number of modes not set." << endl;
    return;
  }

  if (!recalc_needed())
    return;

  if (global.polarisation == TE_TM)
  {
    const int n = int(global.N/2);

    if (2*n != global.N)
      cout << "Warning: changing N to even number." << endl;

    const Complex old_beta = global_slab.beta;

    global.N = n;
    global_slab.beta = 0.0;
    global.polarisation = TE;
    find_modes();

    vector<Mode*> TE_modeset;
    for (unsigned int i=0; i<modeset.size(); i++)
    {
      Slab_M_Mode* mode = new Slab_M_Mode(TE, modeset[i]->get_kz(), this);
      mode->normalise();
      TE_modeset.push_back(mode);
    }

    last_lambda = 0.0; // Force a recalc.

    global.polarisation = TM;
    find_modes();

    modeset.insert(modeset.begin(), TE_modeset.begin(), TE_modeset.end());

    global.N = 2*n;
    global_slab.beta = old_beta;
    global.polarisation = TE_TM;

    return;
  }

  // If we already calculated modes for a different wavelength/gain
  // combination, use these as an initial estimate, else find them
  // from scratch.

  if (global.sweep_from_previous && (modeset.size() == global.N))
    find_modes_by_sweep();
  else
    if (global.solver == ADR)
      find_modes_from_scratch_by_ADR();
    else
      find_modes_from_scratch_by_track();

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_modes_from_scratch_by_ADR
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes_from_scratch_by_ADR()
{
  // Determine reflection coefficients of walls.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  Complex R_left  = l_wall ? l_wall->get_R12() : -1;
  Complex R_right = r_wall ? r_wall->get_R12() : -1;
  
  // Find and max/min real refractive indices of structure.

  Real max_n = 0.0;
  vector<Real> n_lossless;
  for (unsigned int i=0; i<materials.size(); i++)
  {
    Real n = real( materials[i]->n() * sqrt( materials[i]->mur()) );

    n_lossless.push_back(n);
    
    if (n > max_n)
      max_n = n;
  }

  Real min_n = max_n;
  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Real n = n_lossless[i];
    
    if (n < min_n)
      min_n = n;
  }

  // Locate zeros starting from initial contour.

  const Real k0 = 2*pi/global.lambda;
  
  Real max_kt = abs(k0*sqrt(Complex(max_n*max_n - min_n*min_n)))+1;

  Real Re=real(global.C_upperright);
  Real Im=imag(global.C_upperright);

  Complex lowerleft (-0.1,      -0.1);
  Complex upperright( Re*max_kt, Im*max_kt);

  SlabDisp disp(materials, thicknesses, global.lambda, l_wall, r_wall);
  vector<Complex> kt = N_roots(disp, global.N, lowerleft, upperright);
  
  cout << "Calls to slab dispersion relation : "<< disp.times_called() << endl;
  
  // Eliminate false zeros.

  const Real eps_copies = 1e-6;
  
  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Complex kt_i = k0*sqrt(materials[i]->n()*materials[i]->n() - min_n*min_n);
    
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
  
  // Check if enough modes are found.

  if (kt.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kt.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }
  
  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<kt.size(); i++)
  { 
    Complex kz = sqrt(k0*k0*min_n*min_n - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Slab_M_Mode *newmode = new Slab_M_Mode(global.polarisation, kz, this);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_modes_from_scratch_by_track
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes_from_scratch_by_track()
{
  // Set constants
  
  const Real eps        = 1e-13;
  const Real eps_zero   = 1e-10;
  const Real eps_copies = 1e-6;
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;
  
  // Determine reflection coefficients of walls.

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  Complex R_left  = l_wall ? l_wall->get_R12() : -1;
  Complex R_right = r_wall ? r_wall->get_R12() : -1;

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

  // Find and max/min real refractive indices of structure.

  Real max_n = 0.0;
  vector<Real> n_lossless;
  for (unsigned int i=0; i<materials.size(); i++)
  {
    Real n = real( materials[i]->n() * sqrt( materials[i]->mur()) );

    n_lossless.push_back(n);
    
    if (n > max_n)
      max_n = n;
  }

  Real min_n = max_n;
  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Real n = n_lossless[i];
    
    if (n < min_n)
      min_n = n;
  }
  
  // Create dispersion relation for lossless structure.
  
  SlabDisp disp(materials, thicknesses, lambda, l_wall, r_wall);
  params = disp.get_params();
  
  SlabDisp disp_lossless(materials, thicknesses, lambda, l_wall, r_wall);
  vector<Complex> params_lossless;
  for (unsigned int i=0; i<params.size(); i++)
    params_lossless.push_back(real(params[i]));
  disp_lossless.set_params(params_lossless);

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
  
  Wrap_imag_to_abs prop_wrap(disp_lossless);

  Real prop_kt_end_lossless
    = abs(k0*sqrt(Complex(max_n*max_n - min_n*min_n)));

  vector<Real> kt_prop_lossless;
  if (abs(prop_kt_end_lossless) > 0)
  {
    if (d_kt > prop_kt_end_lossless)
      d_kt = prop_kt_end_lossless;

    int sec = (global.precision_enhancement == 1) ? 1 : 0; // security level.

    kt_prop_lossless = brent_all_minima
      (prop_wrap,0.0001,prop_kt_end_lossless,d_kt/global.precision,eps,sec);

    // If precision_enhancement is larger than 1, refine our search for
    // the guided modes.

    if (global.precision_enhancement > 1)
    {
      vector<Real> kt_fine = brent_refine_minima
        (prop_wrap,kt_prop_lossless,global.dx_enhanced,
         d_kt/global.precision/global.precision_enhancement,eps,1);
      kt_prop_lossless = kt_fine;
    }
  }
  else
    prop_kt_end_lossless = 0.25; // To make rectangle for complex zero search.
  
  reverse(kt_prop_lossless.begin(), kt_prop_lossless.end());

  // Eliminate false zeros.

  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Real kt_i = k0*sqrt(n_lossless[i]*n_lossless[i] - min_n*min_n);

    remove_elems(&kt_prop_lossless,  kt_i, eps_copies);
    remove_elems(&kt_prop_lossless, -kt_i, eps_copies);   
    remove_elems(&kt_prop_lossless,  0.0,  eps_copies);
  }
  
  // Check if the minima found correspond really to zeros and increase
  // precision using a mueller solver (in the absense of branchcuts).
  // Don't do this if precision_enhancement is larger than 1, since this
  // typically indicate the presence of nearly generate modes.

  vector<Complex> kt_lossless;
  for (unsigned int i=0; i<kt_prop_lossless.size(); i++)
  { 
    const Real fx = abs(disp_lossless(I*kt_prop_lossless[i]));
    
    if (fx > 1e-2)
    {
      cout << "Warning: possibly insufficient precision around kt "
           << kt_prop_lossless[i] << "." << endl;
      cout << "Removing this candidate with abs(fx) " << fx << "." << endl;
    }
    else
    {
      Real kt_new;
      
      if ( (!branchcut) && (global.precision_enhancement == 1) )
        kt_new = imag(mueller(disp_lossless, I*kt_prop_lossless[i],
                              I*kt_prop_lossless[i]+0.002));
      else
        kt_new = kt_prop_lossless[i];

      kt_lossless.push_back(I*kt_new);
    }
  }


  
  //
  // II: Find evanescent modes of lossless structure.
  //

  if (!(branchcut || TBC))
  {
    Wrap_real_to_abs evan_wrap(disp_lossless);
    
    Real kt_begin = .0001;
    int extra = 1;
    int modes_left = modes_left = (kt_lossless.size() >= global.N + 2) ? 0
      : global.N - kt_lossless.size() + 2;
    
    while (modes_left)
    {
      if (extra > 1)
        cout << "Lost too many candidate modes. Expanding search region."
             << endl;
      
      vector<Real> kt_evan_lossless = brent_N_minima
        (evan_wrap,kt_begin,modes_left+extra,d_kt/global.precision_rad,eps,1);
      
      // Check if the minima found correspond really to zeros and increase
      // precision with mueller solver.

      for (unsigned int i=0; i<kt_evan_lossless.size(); i++)
      { 
        const Real fx = abs(disp_lossless(kt_evan_lossless[i]));
        if (fx > 1e-2)
        {
          cout << "Warning: possibly insufficient precision around kt "
               <<  kt_evan_lossless[i] << "." << endl;
          cout << "Removing this candidate with abs(fx) "
               << fx << "." << endl;
        }
        else
        {
          Real kt_new = branchcut
            ? kt_evan_lossless[i]
            : real(mueller(disp_lossless,kt_evan_lossless[i],
                           kt_evan_lossless[i]+0.002*I));

          kt_lossless.push_back(kt_new);
        }
      }

      modes_left = (kt_lossless.size() >= global.N + 2) ? 0
        : global.N - kt_lossless.size() + 2;

      kt_begin = real(kt_lossless[kt_lossless.size()-1])+.0001;
      extra++;
    }
  }

  
  
  //
  // III: Eliminate doubles and flase zeros and trace modes to those of
  //      true structure.
  //

  if (global.precision_enhancement == 1)
    remove_copies(&kt_lossless, eps_copies);

  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Complex kt_i = k0*sqrt(n_lossless[i]*n_lossless[i] - min_n*min_n);

    remove_elems(&kt_lossless,     kt_i,     eps_copies);
    remove_elems(&kt_lossless,    -kt_i,     eps_copies);   
    remove_elems(&kt_lossless, Complex(0.0), eps_copies);
  }

  vector<Complex> kt, forbidden;
  if (global.chunk_tracing == true)
    kt = traceroot_chunks
      (kt_lossless, disp_lossless, disp, forbidden, global.sweep_steps);
  else
    kt = traceroot
      (kt_lossless, disp_lossless, disp, forbidden, global.sweep_steps);

 

  //
  // IV: Find modes in complex plane
  //

  if ( (branchcut || TBC) && (kt.size() < global.N) )
  {
    Real evan_kt_end = 2*(global.N-kt.size())*d_kt; // Not general!

    Real Re=real(global.C_upperright);
    Real Im=imag(global.C_upperright);

    Complex lowerleft (0.001,          0.001);
    //Complex upperright(Re*evan_kt_end, Im*prop_kt_end_lossless);
    Complex upperright = global.C_upperright;

    const int sections = int(global.C_steps);

    vector<Complex> kt_complex = allroots(disp, lowerleft, upperright);

    // Eliminate doubles and false zeros.

    remove_copies(&kt_complex, eps_copies);

    for (unsigned int i=0; i<n_lossless.size(); i++)
    {
      Complex kt_i = k0*sqrt(n_lossless[i]*n_lossless[i] - min_n*min_n);

      remove_elems(&kt_complex,  kt_i, eps_copies);
      remove_elems(&kt_complex, -kt_i, eps_copies);
    }

    cout << "Found " << kt_complex.size() << " complex zeros in region "
         << lowerleft << " " << upperright << "." <<endl;

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
  
  // Check if enough modes are found.

  if (kt.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kt.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }
  
  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<kt.size(); i++)
  { 
    Complex kz = sqrt(k0*k0*min_n*min_n - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Slab_M_Mode *newmode = new Slab_M_Mode(global.polarisation, kz, this);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab_M::find_modes_by_sweep
//
/////////////////////////////////////////////////////////////////////////////

void Slab_M::find_modes_by_sweep()
{
  // Set constants.
  
  const Real k0 = 2*pi/global.lambda;

  SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;
    
  // Trace modes from old configuration to new one.
  
  SlabDisp disp_new(materials, thicknesses, global.lambda, l_wall, r_wall);
  SlabDisp disp_old(materials, thicknesses, global.lambda, l_wall, r_wall);
  disp_old.set_params(params);

  Complex min_n = disp_old.get_min_n();
  Complex k0_old = 2*pi/last_lambda;
  
  vector<Complex> kt_old;
  for (unsigned int i=0; i<modeset.size(); i++)
    kt_old.push_back(sqrt(pow(k0_old*min_n, 2)-pow(modeset[i]->get_kz(), 2)));
  
  vector<Complex> forbidden;
  vector<Complex> kt =
    traceroot(kt_old, disp_old, disp_new, forbidden, global.sweep_steps);

  params = disp_new.get_params();

  // Check if modes were lost during tracing.
  
  if (kt.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kt.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }

  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<kt.size(); i++)
  {
    Complex kz = sqrt(k0*k0*min_n*min_n - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Slab_M_Mode *newmode = new Slab_M_Mode(global.polarisation, kz, this);
      
    newmode->normalise();

    modeset.push_back(newmode);
  }
  
  sort_modes();
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
    if ( real(materials[i]->n()) > real(materials[core_index]->n()) )
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
  if (global.polarisation == TE_TM)
  {
    const int n = int(global.N/2);

    if (2*n != global.N)
      cout << "Warning: changing N to even number." << endl;

    const Complex old_beta = global_slab.beta;

    global.N = n;
    global_slab.beta = 0.0;
    global.polarisation = TE;
    find_modes_single_pol();

    vector<Mode*> TE_modeset;
    for (unsigned int i=0; i<modeset.size(); i++)
    {
      UniformSlabMode* mode 
        = new UniformSlabMode(TE, modeset[i]->get_kz(), this);
      mode->normalise();
      TE_modeset.push_back(mode);
    }

    last_lambda = 0.0; // Force a recalc.

    global.polarisation = TM;
    find_modes_single_pol();

    modeset.insert(modeset.begin(), TE_modeset.begin(), TE_modeset.end());

    global.N = 2*n;
    global_slab.beta = old_beta;
    global.polarisation = TE_TM;

    return;
  }

  find_modes_single_pol();
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformSlab::find_modes_single_pol
//
/////////////////////////////////////////////////////////////////////////////

void UniformSlab::find_modes_single_pol()
{
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
    cout << "Error: wavelength not set." << endl;
    return;
  }
  
  if (global.N == 0)
  {
    cout << "Error: number of modes not set." << endl;
    return;
  }

  // Set contants and clear modeset.
  
  const Complex k = 2*pi/global.lambda * core->n();

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // If no analytic solution is available, follow same route as
  // non-uniform slab.

  if (!analytic)
  {
    Slab_M tmp = (*core)(get_width()/2.) + (*core)(get_width()/2.);

    if (l_wall)
      tmp.set_left_wall(*l_wall);

    if (r_wall)
      tmp.set_right_wall(*r_wall);
    
    tmp.find_modes();

    for (unsigned int i=1; i<=tmp.N(); i++)
    {
      UniformSlabMode *newmode = new 
        UniformSlabMode(global.polarisation, tmp.get_mode(i)->get_kz(), this);

      newmode->normalise();

      modeset.push_back(newmode);
    }

    // Remember wavelength and gain these modes were calculated for.

    last_lambda = global.lambda;

    if (global.gain_mat)
      last_gain_mat = *global.gain_mat;

    return;
  }

  // Use analytical solution available. 

  // TEM mode.
  
  if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
                                     && (abs(R_right + 1.0) < 1e-10) )
    || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
                                     && (abs(R_right - 1.0) < 1e-10) ) )
  {
    UniformSlabMode *newmode
      = new UniformSlabMode(global.polarisation, k, this);

    newmode->normalise();
    
    modeset.push_back(newmode);
  }

  // Other modes.
  
  Complex start  = (abs(R_left*R_right-1.0) < 1e-10) ? 0 : pi;
  unsigned int i = (abs(R_left*R_right-1.0) < 1e-10) ? 1 : 0;
  
  while (modeset.size() != global.N)
  {
    Complex kz = sqrt(k*k - pow((start+2*i*pi)/2./get_width(), 2));

    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    UniformSlabMode *newmode
      = new UniformSlabMode(global.polarisation, kz, this);

    newmode->normalise();
    
    modeset.push_back(newmode);

    i++;
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
