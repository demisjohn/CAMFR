
/////////////////////////////////////////////////////////////////////////////
//
// File:     circ.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990824
// Version:  1.2
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include "circ.h"
#include "circdisp.h"
#include "circmode.h"
#include "circoverlap.h"
#include "../../math/calculus/calculus.h"
#include "../../math/calculus/quadrature/patterson_quad.h"
#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// CircGlobal
//
/////////////////////////////////////////////////////////////////////////////

CircGlobal global_circ = {1, cos_type};



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::hankel
//
/////////////////////////////////////////////////////////////////////////////

Hankel Circ_M::hankel = kind_1;


/////////////////////////////////////////////////////////////////////////////
//
// CircFlux
//
/////////////////////////////////////////////////////////////////////////////

class CircFlux : public RealFunction
{
  public:

    CircFlux(const FieldExpansion& fe_) : fe(fe_) {}

    Real operator()(const Real& rho)
    {
      counter++;
      Field f=fe.field(Coord(rho,0,0));
      return real( pi*rho*(f.E1*conj(f.H2)-f.E2*conj(f.H1)) );
    }

  protected:

    FieldExpansion fe;
};



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::S_flux()
//
/////////////////////////////////////////////////////////////////////////////

Real Circ_M::S_flux(const FieldExpansion& f,
                    Real c1_start, Real c1_stop,
                    Real precision = 1e-10) const
{  
  CircFlux flux(f);
  return patterson_quad(flux, c1_start, c1_stop, precision);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::createStretcher
//
/////////////////////////////////////////////////////////////////////////////

CoordStretcher Circ_M::createStretcher()
{ 
  Real c_start = real(radius[M-2] + radius[M-1]) / 2; // halfway outer ring
  
  return CoordStretcher(c_start,   (     radius[M-1]  - c_start)
                                 / (real(radius[M-1]) - c_start));
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::operator==
//
/////////////////////////////////////////////////////////////////////////////

bool Circ_M::operator==(const Waveguide& w_) const
{
  const Circ_M* w = dynamic_cast<const Circ_M*>(&w_);
  
  for (unsigned int i=0; i<material.size(); i++)
    if (material[i] != w->material[i])
      return false;

  for (unsigned int i=0; i<radius.size(); i++)
    if ( abs(radius[i] - w->radius[i]) > 1e-12)
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::contains
//
/////////////////////////////////////////////////////////////////////////////

bool Circ_M::contains(const Material& m) const
{
  for (unsigned int i=0; i<material.size(); i++)
    if (material[i] == &m)
      return true;
  
  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool Circ_M::no_gain_present() const
{
  for (unsigned int i=0; i<material.size(); i++)
    if (!material[i]->no_gain_present())
      return false;
  
  return true;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::find_modes()
//
/////////////////////////////////////////////////////////////////////////////

void Circ_M::find_modes()
{
  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  cerr << "Circ_M::findModes to be implemented." << endl;
  exit (-1);
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_M::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

void Circ_M::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL)
{  
  const Circ_M* medium_I  = this;
  const Circ_M* medium_II = dynamic_cast<const Circ_M*>(w);

  // Make sorted list of evaluation points for field cache.
  
  vector<Complex> disc(medium_I->radius);

  disc.push_back(0.0);

  for (unsigned int k=0; k<medium_II->radius.size(); k++)
    disc.push_back(medium_II->radius[k]);

  remove_copies(&disc, 1e-6);
    
  sort(disc.begin(), disc.end(), RealSorter());

  // Fill field cache.

  CircCache cache(global.N, disc.size()-1);

  for (int i=1; i<=int(global.N); i++)
  {
    const Circ_M_Mode* mode_I
      = dynamic_cast<const Circ_M_Mode*>(medium_I ->get_mode(i));

    const Circ_M_Mode* mode_II
      = dynamic_cast<const Circ_M_Mode*>(medium_II->get_mode(i));

    for (int k=0; k<int(disc.size()-1); k++)
    {
      const bool ang_dep = false;
      
      const Coord lower(disc[k],  0,0, Plus);
      const Coord upper(disc[k+1],0,0, Min);

      Complex dE_I_l, dH_I_l, dE_II_l, dH_II_l;
      Complex dE_I_u, dH_I_u, dE_II_u, dH_II_u;

      Field f_I_l  = mode_I ->field_at(lower, &dE_I_l,  &dH_I_l,  ang_dep);
      Field f_II_l = mode_II->field_at(lower, &dE_II_l, &dH_II_l, ang_dep);
      Field f_I_u  = mode_I ->field_at(upper, &dE_I_u,  &dH_I_u,  ang_dep);
      Field f_II_u = mode_II->field_at(upper, &dE_II_u, &dH_II_u, ang_dep);

      cache. E_l(1,i,k+1) = f_I_l.Ez;
      cache. H_l(1,i,k+1) = f_I_l.Hz;
      cache. E_u(1,i,k+1) = f_I_u.Ez;
      cache. H_u(1,i,k+1) = f_I_u.Hz;
      cache.dE_l(1,i,k+1) = dE_I_l;
      cache.dH_l(1,i,k+1) = dH_I_l;
      cache.dE_u(1,i,k+1) = dE_I_u;
      cache.dH_u(1,i,k+1) = dH_I_u;

      cache. E_l(2,i,k+1) = f_II_l.Ez;
      cache. H_l(2,i,k+1) = f_II_l.Hz;
      cache. E_u(2,i,k+1) = f_II_u.Ez;
      cache. H_u(2,i,k+1) = f_II_u.Hz;
      cache.dE_l(2,i,k+1) = dE_II_l;
      cache.dH_l(2,i,k+1) = dH_II_l;
      cache.dE_u(2,i,k+1) = dE_II_u;
      cache.dH_u(2,i,k+1) = dH_II_u;
    }
  }

  // Calculate overlap matrices.

  for (int i=1; i<=int(global.N); i++)
    for (int j=1; j<=int(global.N); j++)
    {      
      (*O_I_II)(i,j) = overlap
        (dynamic_cast<const Circ_M_Mode*>(medium_I ->get_mode(i)),
         dynamic_cast<const Circ_M_Mode*>(medium_II->get_mode(j)),
         &cache, &disc, i, j, 1, 2);
      
      (*O_II_I)(i,j) = overlap
        (dynamic_cast<const Circ_M_Mode*>(medium_II->get_mode(i)),
         dynamic_cast<const Circ_M_Mode*>(medium_I ->get_mode(j)),
         &cache, &disc, i, j, 2, 1);

      if (O_I_I) (*O_I_I)(i,j) = overlap
        (dynamic_cast<const Circ_M_Mode*>(medium_I ->get_mode(i)),
         dynamic_cast<const Circ_M_Mode*>(medium_I ->get_mode(j)),
         &cache, &disc, i, j, 1, 1);

      if (O_II_II) (*O_II_II)(i,j) = overlap
        (dynamic_cast<const Circ_M_Mode*>(medium_II->get_mode(i)),
         dynamic_cast<const Circ_M_Mode*>(medium_II->get_mode(j)),
         &cache, &disc, i, j, 2, 2);
    }
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2::Circ_2
//
/////////////////////////////////////////////////////////////////////////////

Circ_2::Circ_2(const Complex& r, Material& core_,
               const Complex& R, Material& cladding_)
{
  radius.push_back(r);
  radius.push_back(R);
  material.push_back(&core_);
  material.push_back(&cladding_);
  M       = 2;
  uniform = false;
  core    = &core_;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2::find_modes()
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

  const bool lossy =
       imag(radius[0])          || imag(radius[1])
    || imag(material[0]->n())   || imag(material[1]->n())
    || imag(material[0]->mur()) || imag(material[1]->mur());

  // If we already calculated modes for a different wavelength/gain
  // combination, use these as an initial estimate, else find them
  // from scratch.
  
  if (global.sweep_from_previous && lossy && (modeset.size() == global.N))
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
// Circ_2::find_modes_from_scratch_by_ADR
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2::find_modes_from_scratch_by_ADR()
{
  // Set constants.
  
  const Real k0         = 2*pi/global.lambda;

  const Complex   n1    = material[0]->n();
  const Complex   n2    = material[1]->n();
  const Complex  mur1   = material[0]->mur();
  const Complex  mur2   = material[1]->mur(); 
  const Complex delta_k = n2*n2*mur2 - n1*n1*mur1; // neg if guiding
  
  const Complex   k1    = k0*n1*sqrt(mur1);
  const Complex   k2    = k0*n2*sqrt(mur2);
  
  // Locate zeros starting from initial contour.

  Real max_kt = abs(k0*sqrt(Complex(
    real(n1)*real(n1)*real(mur1) - real(n2)*real(n2)*real(mur2) )));

  Real Re=real(global.C_upperright);
  Real Im=imag(global.C_upperright);

  Complex lowerleft (-Re*max_kt,   -0.1);
  Complex upperright( Re*max_kt, Im*max_kt);

  const bool scale_complex = true;
  
  Circ_2_closed disp
    (radius[0], radius[1], *material[0], *material[1], global.lambda,
     global_circ.order, rad, hankel, global.polarisation, scale_complex);
  
  vector<Complex> kr2 = N_roots(disp, global.N, lowerleft, upperright);
  
  cout << "Calls to circ dispersion relation : "<< disp.times_called() << endl;
  
  // Eliminate copies and false zeros.

  const Real eps_copies = 1e-6;

  remove_copies(&kr2, eps_copies);
  
  remove_elems(&kr2,   Complex(0.0),     eps_copies);  
  remove_elems(&kr2, + k0*sqrt(delta_k), eps_copies);
  remove_elems(&kr2, - k0*sqrt(delta_k), eps_copies);
  
  // Check if enough modes are found.

  if (kr2.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kr2.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }
  
  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  for (unsigned int i=0; i<kr2.size(); i++)
  {
    Complex kz  = signedsqrt(k2*k2 - kr2[i]*kr2[i], material[0]);
    Complex kr1 = signedsqrt(k1*k1 -     kz*kz,     material[0]);

    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    if (real(kr2[i]) < -.001) // Backward mode of complex set.
      kz = -kz;
    
    Circ_2_Mode *newmode
      = new Circ_2_Mode(Polarisation(unknown), kz, kr1, kr2[i], this);

    newmode->normalise();

    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2::find_modes_from_scratch_by_track
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2::find_modes_from_scratch_by_track()
{
  // Set constants.
  
  const Real eps        = 1e-13;
  const Real eps_copies = 1e-6;
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;

  const Complex   n1    = material[0]->n();
  const Complex   n2    = material[1]->n();
  const Complex  mur1   = material[0]->mur();
  const Complex  mur2   = material[1]->mur(); 
  const Complex delta_k = n2*n2*mur2 - n1*n1*mur1; // neg if guiding
  
  const Complex   k1    = k0*n1*sqrt(mur1);
  const Complex   k2    = k0*n2*sqrt(mur2);

  //
  // (I) : Find kr2 of guided modes (along pos. imag. kr2 axis).
  //

  // Find modes for lossless structure.

  Real guided_dk = abs(0.5*pi*lambda/real(radius[0])); // roughly
  Real guided_k_end_lossless = abs(k0*sqrt(Complex(
    real(n1)*real(n1)*real(mur1) - real(n2)*real(n2)*real(mur2) )));

  Material co_lossless(real(n1), real(mur1));
  Material cl_lossless(real(n2), real(mur2));
  
  Circ_2_closed guided_disp_lossless
    (radius[0], real(radius[1]), co_lossless, cl_lossless, lambda,
     global_circ.order, guided, hankel, global.polarisation);

  Wrap_imag_to_real guided_wrap(guided_disp_lossless);
    
  vector<Real> kr2_guided_lossless_real;

  if (abs(delta_k) > 0) // else homogeneously filled and no modes here
  {
    if (guided_dk > guided_k_end_lossless)
      guided_dk = guided_k_end_lossless;
    
    kr2_guided_lossless_real
      = brent_all_roots(guided_wrap,0.0001,guided_k_end_lossless,
                        guided_dk/global.precision,eps,1);
  }

  reverse(kr2_guided_lossless_real.begin(), kr2_guided_lossless_real.end());

  vector<Complex> kr2_guided_lossless;
  for (unsigned int i=0; i<kr2_guided_lossless_real.size(); i++)
    kr2_guided_lossless.push_back(I*kr2_guided_lossless_real[i]);
  
  // Trace modes to those of true structure.

  Circ_2_closed guided_disp
    (radius[0], radius[1], *material[0], *material[1],
     lambda, global_circ.order, guided, hankel, global.polarisation);

  vector<Complex> forbidden;
  forbidden.push_back(0.0);
  forbidden.push_back(+k0*sqrt(delta_k));
  forbidden.push_back(-k0*sqrt(delta_k));
  
  vector<Complex> kr2_guided =
    traceroot(kr2_guided_lossless,guided_disp_lossless,guided_disp,
              forbidden, global.sweep_steps);

  guided_disp_params = guided_disp.get_params();
  
  cout << "Calls to guided dispersion relation: "
       << guided_disp_lossless.times_called() + guided_disp.times_called()
       << endl;


  
  //
  // (II) : Find kr2 of radiation modes (along pos. real kr2 axis).
  //

  // Find modes for lossless structure.

  Real rad_dk = (global_circ.order == 0) ? pi/real(radius[0])
                                         : 0.5*pi/real(radius[0]); // roughly

  Circ_2_closed_rad_lossless rad_disp_lossless
    (radius[0], real(radius[1]), co_lossless, cl_lossless,
     lambda, global_circ.order, hankel, global.polarisation);
  
  Wrap_real_to_real rad_wrap(rad_disp_lossless);
  
  int rad_modes = global.N - kr2_guided.size() + 2; // false zeros

  vector<Real> kr2_rad_lossless;
  Real rad_k_end_lossless;
  
  if (rad_modes <= 0)
  {
    cout << "Warning: no radiation modes included." << endl;
    rad_k_end_lossless = 0.0;
  }
  else
  {
    kr2_rad_lossless = brent_N_roots
      (rad_wrap, 0.0001, rad_modes, rad_dk/global.precision_rad, eps, 1);
    rad_k_end_lossless = kr2_rad_lossless[kr2_rad_lossless.size()-1];
  }
  
  // Trace modes to those of true structure.
  
  Circ_2_closed rad_disp
    (radius[0], radius[1], *material[0], *material[1],
     lambda, global_circ.order, rad, hankel, global.polarisation);

  for (unsigned int i=0; i<kr2_guided.size(); i++)
    forbidden.push_back(kr2_guided[i]);
  
  vector<Complex> kr2_rad;

  if (global.chunk_tracing == true)
    kr2_rad = traceroot_chunks
      (kr2_rad_lossless, rad_disp_lossless, rad_disp,
       forbidden, global.sweep_steps);
  else
    kr2_rad = traceroot
      (kr2_rad_lossless, rad_disp_lossless, rad_disp,
       forbidden, global.sweep_steps);

  for (unsigned int i=0; i<kr2_rad.size(); i++)
    if (    ( (imag(kr2_rad[i]) < 0) && (hankel == kind_1) )
         || ( (imag(kr2_rad[i]) > 0) && (hankel == kind_2) ) )
      kr2_rad[i] = -kr2_rad[i]; // Pick stable sign.
  
  rad_disp_params = rad_disp.get_params();
  
  cout << "Calls to radiation dispersion relation: "
       << rad_disp_lossless.times_called() + rad_disp.times_called()
       << endl;

  
/*  
  //
  // (III) : Find kr2 of complex modes (rare modes not lying on the axes).
  //

  // Find modes for lossless structure.
    
  const bool scale_complex = true;
  
  Circ_2_closed_cutoff complex_disp_lossless
    (radius[0], real(radius[1]), co_lossless, cl_lossless, lambda,
     global_circ.order, rad, hankel, global.polarisation, scale_complex);

  if (abs(guided_k_end_lossless) == 0)
    guided_k_end_lossless = 1; // To make a rectangle.

  Real Re=real(global.C_upperright);
  Real Im=imag(global.C_upperright);

  Complex lowerleft (0.001,                 0.001);
  Complex upperright(Re*rad_k_end_lossless, Im*guided_k_end_lossless);
  
  vector<Complex> kr2_complex_lossless_try = allroots
    (complex_disp_lossless, lowerleft, upperright);
  
  vector<Complex> kr2_complex_lossless;
  for (unsigned int i=0; i<kr2_complex_lossless_try.size(); i++)
    if (    ( abs(real(kr2_complex_lossless_try[i])) < eps_copies)
         || ( abs(imag(kr2_complex_lossless_try[i])) < eps_copies) )
    cout << "Warning: unable to locate possible complex mode." << endl;
  else
    kr2_complex_lossless.push_back(kr2_complex_lossless_try[i]);
  
  if (kr2_complex_lossless.size() > 0)
    cout << "Found " << kr2_complex_lossless.size()
         << " set(s) of complex modes!" << endl;

  // Add other complex mode of symmetric quartet.
  
  vector<Complex> kr2_complex_lossless_bis = kr2_complex_lossless;
  for (unsigned int i=0; i<kr2_complex_lossless_bis.size(); i++)
    kr2_complex_lossless.push_back
      (Complex(-real(kr2_complex_lossless_bis[i]),
               +imag(kr2_complex_lossless_bis[i])));
  
  // Trace modes to those of true structure.
  
  Circ_2_closed_cutoff complex_disp
    (radius[0], radius[1], *material[0], *material[1],lambda,
     global_circ.order, rad, hankel, global.polarisation, scale_complex);

  for (unsigned int i=0; i<kr2_rad.size(); i++)
    forbidden.push_back(kr2_rad[i]);
  
  vector<Complex> kr2_complex = traceroot
    (kr2_complex_lossless, complex_disp_lossless, complex_disp,
     forbidden, global.sweep_steps);

  // Remember second half of complex modes as backward modes.
  
  for (unsigned int i=int(kr2_complex.size()/2); i<kr2_complex.size(); i++)
    kr2_backward.push_back(kr2_complex[i]);
 
  cout << "Calls to complex dispersion relation: "
       << complex_disp_lossless.times_called() + complex_disp.times_called()
       << endl;
*/ 
  //
  // Join modes, eliminate double zeros and false zeros at kr2=0 and kr1=0.
  //
  
  vector<Complex> kr2 = kr2_guided;

  //for (unsigned int i=0; i<kr2_complex.size(); i++)
  //  kr2.push_back(kr2_complex[i]);
  
  for (unsigned int i=0; i<kr2_rad.size(); i++)
    kr2.push_back(kr2_rad[i]);
  
  remove_copies(&kr2, eps_copies);
  
  remove_elems(&kr2,   Complex(0.0),     eps_copies);  
  remove_elems(&kr2, + k0*sqrt(delta_k), eps_copies);
  remove_elems(&kr2, - k0*sqrt(delta_k), eps_copies);

  // Check if enough modes are found. If needed, eliminate the highest
  // order radiation modes.
  
  if (kr2.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kr2.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }
  else
    kr2.erase(kr2.begin() + global.N, kr2.end());

  
  
  //
  // Create modeset.
  //

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  no_of_guided_modes = 0;
  
  for (unsigned int i=0; i<kr2.size(); i++)
  {
    Complex kz  = signedsqrt(k2*k2 - kr2[i]*kr2[i], material[0]);
    Complex kr1 = signedsqrt(k1*k1 -     kz*kz,     material[0]);

    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    //if (find(kr2_backward.begin(), kr2_backward.end(), kr2[i])
    //    != kr2_backward.end())
    //  kz = -kz;
    
    Circ_2_Mode *newmode
      = new Circ_2_Mode(Polarisation(unknown), kz, kr1, kr2[i], this);

    newmode->normalise();

    modeset.push_back(newmode);

    if (abs(real(kz)) > abs(imag(kz)))
      no_of_guided_modes++;
  }
  
  sort_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_2::find_modes_by_sweep
//
/////////////////////////////////////////////////////////////////////////////

void Circ_2::find_modes_by_sweep()
{ 
  // Set constants.
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;

  const Complex   n1    = material[0]->n();
  const Complex   n2    = material[1]->n();
  const Complex  mur1   = material[0]->mur();
  const Complex  mur2   = material[1]->mur(); 
  const Complex delta_k = n2*n2*mur2 - n1*n1*mur1; // neg if guiding
  
  const Complex   k1    = k0*n1*sqrt(mur1);
  const Complex   k2    = k0*n2*sqrt(mur2);

  // Guided modes.

  vector<Complex> kr2_guided_old;
  for (unsigned int i=0; i<no_of_guided_modes; i++)
  {
    const Circ_M_Mode* mode = dynamic_cast<const Circ_M_Mode*>(modeset[i]);
    kr2_guided_old.push_back(mode->kr_at(radius[1]));
  }

  Circ_2_closed
    guided_disp_new(radius[0], radius[1], *material[0], *material[1],
                    lambda, global_circ.order, guided, hankel,
                    global.polarisation);

  Circ_2_closed
    guided_disp_old(radius[0], radius[1], *material[0], *material[1],
                    lambda, global_circ.order, guided, hankel,
                    global.polarisation);

  guided_disp_old.set_params(guided_disp_params);

  vector<Complex> forbidden;
  forbidden.push_back(0.0);
  forbidden.push_back(+k0*sqrt(delta_k));
  forbidden.push_back(-k0*sqrt(delta_k));

  vector<Complex> kr2_guided =
    traceroot(kr2_guided_old, guided_disp_old, guided_disp_new,
              forbidden, global.sweep_steps);

  guided_disp_params = guided_disp_new.get_params();

  // Radiation and complex modes.

  vector<Complex> kr2_rad_old;
  for (unsigned int i=no_of_guided_modes; i<modeset.size(); i++)
  {
    const Circ_M_Mode* mode = dynamic_cast<const Circ_M_Mode*>(modeset[i]);
    kr2_rad_old.push_back(mode->kr_at(radius[1]));
  }

  Circ_2_closed
    rad_disp_new(radius[0], radius[1], *material[0], *material[1],
                 lambda, global_circ.order, rad, hankel, global.polarisation);

  Circ_2_closed
    rad_disp_old(radius[0], radius[1], *material[0], *material[1],
                 lambda, global_circ.order, rad, hankel, global.polarisation);

  rad_disp_old.set_params(rad_disp_params);

  for (unsigned int i=0; i<kr2_guided.size(); i++)
    forbidden.push_back(kr2_guided[i]);

  vector<Complex> kr2_rad = 
    traceroot(kr2_rad_old, rad_disp_old, rad_disp_new,
              forbidden, global.sweep_steps);

  rad_disp_params = rad_disp_new.get_params();

  // Backward modes. These are also included in the radiation modes,
  // but we need to keep track of them seperately.

  vector<Complex> kr2_backward_new =
      traceroot(kr2_backward, rad_disp_old, rad_disp_new,
                forbidden, global.sweep_steps);
  kr2_backward = kr2_backward_new;

  // Join modes.

  vector<Complex> kr2 = kr2_guided;

  for (unsigned int i=0; i<kr2_rad.size(); i++)
    kr2.push_back(kr2_rad[i]);

  // Check if modes were lost during tracing.
  
  if (kr2.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << kr2.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }

  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  for (unsigned int i=0; i<kr2.size(); i++)
  {
    Complex kz  = signedsqrt(k2*k2 - kr2[i]*kr2[i], material[0]);
    Complex kr1 = signedsqrt(k1*k1 -     kz*kz,     material[0]);

    if (real(kz) < 0)
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    if (find(kr2_backward.begin(), kr2_backward.end(), kr2[i])
        != kr2_backward.end())
      kz = -kz;

    Circ_2_Mode *newmode
      = new Circ_2_Mode(Polarisation(unknown), kz, kr1, kr2[i], this);

    newmode->normalise();

    modeset.push_back(newmode);
  }

  sort_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_1::Circ_1
//
/////////////////////////////////////////////////////////////////////////////

Circ_1::Circ_1(const Complex& r, Material& m)
{
  radius.push_back(r);
  material.push_back(&m);
  M       = 1;
  uniform = true;
  core    = &m;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ_1::find_modes()
//
/////////////////////////////////////////////////////////////////////////////

void Circ_1::find_modes()
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
  
  // Set constants.
  
  const Real eps             = 1e-14;
  const Real eps2            = 1e-10;
  const unsigned int N       = global.N;
  const unsigned int N_TE_TM = (global_circ.order != 0) ? int(N/2)+1 : N+1;
  
  // Find TE modes.
  
  Circ_1_TE TE_disp(global_circ.order);
  vector<Real> kr_rho_TE;
  if (     (global_circ.order != 0)
       || ((global_circ.order == 0) && global.polarisation == TE) )
    kr_rho_TE = brent_N_roots(TE_disp, 0, N_TE_TM, pi/2, eps);

  // Find TM modes.

  Circ_1_TM TM_disp(global_circ.order);
  vector<Real> kr_rho_TM;
  if (     (global_circ.order != 0)
       || ((global_circ.order == 0) && global.polarisation == TM) )
    kr_rho_TM = brent_N_roots(TM_disp, 0, N_TE_TM, pi/2, eps);
  
  // Divide by R and join.

  vector<Complex> kr;
  
  for (unsigned int i=0; i<kr_rho_TE.size(); i++)
    if (abs(kr_rho_TE[i]) > eps2) // false zero
      kr.push_back(kr_rho_TE[i] / radius[0]);

  for (unsigned int i=0; i<kr_rho_TM.size(); i++)
    if (abs(kr_rho_TM[i]) > eps2) // false zero
      kr.push_back(kr_rho_TM[i] / radius[0]);

  // Sort and trim.

  sort(kr.begin(), kr.end(), betasorter());
  reverse(kr.begin(), kr.end());
  
  if (kr.size() < N) // Shouldn't happen.
  {
    cerr << "Error: didn't find enough modes ("
         << kr.size() << "/" << global.N << ")."  << endl;
    exit (-1);
  }
  else
    kr.erase(kr.begin() + N, kr.end());
  
  // Create layerset.
  
  const Complex k1_2 = pow(2*pi/global.lambda, 2) * material[0]->epsr()
                                                  * material[0]->mur();

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<kr.size(); i++)
  {
    Complex kz = signedsqrt(k1_2 - kr[i]*kr[i], material[0]);

    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Circ_1_Mode *newmode
      = new Circ_1_Mode(Polarisation(unknown), kz, kr[i], this);

    newmode->normalise();
    
    modeset.push_back(newmode);
  }

  sort_modes(); 

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ::Circ
//
/////////////////////////////////////////////////////////////////////////////

Circ::Circ(const Expression& ex)
{
  // Check expression.
  
  Expression e = ex.flatten();

  for (unsigned int i=0; i<e.get_size(); i++)
    if (!dynamic_cast<Material*>(e.get_term(i)->get_mat()))
    {
      cerr << "Error: expression contains non-material term." << endl;
      return;
    }

  // Build Circ structure.

  if (e.get_size() == 1)
  {
    Material* m = dynamic_cast<Material*>(e.get_term(0)->get_mat());
    Complex   d = e.get_term(0)->get_d();
    
    c = new Circ_1(d, *m);
  }
  else if (e.get_size() == 2)
  {
    Material* m1 = dynamic_cast<Material*>(e.get_term(0)->get_mat());
    Complex   d1 = e.get_term(0)->get_d();

    Material* m2 = dynamic_cast<Material*>(e.get_term(1)->get_mat());
    Complex   d2 = e.get_term(1)->get_d();

    c = new Circ_2(d1, *m1, d1+d2, *m2);
  }
  else
    cout << "Error: more than two materials not yet supported in Circ."
         << endl;

  uniform = c->is_uniform();
  core = c->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Circ::Circ
//
/////////////////////////////////////////////////////////////////////////////

Circ::Circ(const Term& t)
{
  Material* m = dynamic_cast<Material*>(t.get_mat());
  Complex   d = t.get_d();

  if (!m)
  {
    cerr << "Error: expression contains non-material term." << endl;
    return;
  }

  c = new Circ_1(d, *m);

  uniform = c->is_uniform();
  core = c->get_core();
}
