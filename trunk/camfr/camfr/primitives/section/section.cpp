
/////////////////////////////////////////////////////////////////////////////
//
// File:     section.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020129
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "section.h"
#include "sectiondisp.h"
#include "sectionmode.h"
#include "sectionoverlap.h"

using std::vector;
using std::cout;
using std::endl;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// SectionImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

void SectionImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I, cMatrix* O_II_II)
{
  SectionImpl* medium_I  = this;
  SectionImpl* medium_II = dynamic_cast<SectionImpl*>(w);
  
  for (int i=1; i<=int(global.N); i++)
    for (int j=1; j<=int(global.N); j++)
    {
      (*O_I_II)(i,j) = overlap
        (dynamic_cast<const SectionMode*>(medium_I ->get_mode(i)),
         dynamic_cast<const SectionMode*>(medium_II->get_mode(j)));

      (*O_II_I)(i,j) = overlap
        (dynamic_cast<const SectionMode*>(medium_II->get_mode(i)),
         dynamic_cast<const SectionMode*>(medium_I ->get_mode(j)));
      
      if (O_I_I) (*O_I_I)(i,j) = overlap
         (dynamic_cast<const SectionMode*>(medium_I ->get_mode(i)),
          dynamic_cast<const SectionMode*>(medium_I ->get_mode(j)));
      
      if (O_II_II) (*O_II_II)(i,j) = overlap
         (dynamic_cast<const SectionMode*>(medium_II->get_mode(i)),
          dynamic_cast<const SectionMode*>(medium_II->get_mode(j)));
    }
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(const Expression& expression, int M)
{
  // 1D section?

  if (expression.get_size() == 1)
  {
    Slab* slab = dynamic_cast<Slab*>(expression.get_term(0)->get_wg());
    Complex d = expression.get_term(0)->get_d();

    if (slab)
      s = new Section1D(*slab, d);
    else
    {
      py_error("Error: expected a slab to initialise a section.");
      exit (-1);
    }    

    uniform = s->is_uniform();
    core = s->get_core();

    return;
  }
  
  // Transparently factor expression into right and left expression.

  py_error("Automatic splitting of expressions not yet implemented.");
  exit (-1);
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(const Expression& left_ex, const Expression& right_ex, int M)
{
  s = new Section2D(left_ex, right_ex, M);
  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::Section2D
//
/////////////////////////////////////////////////////////////////////////////

Section2D::Section2D(const Expression& left_ex, const Expression& right_ex,
                     int M_)
  : left(left_ex), right(right_ex), M(M_)
{
  // Check values.

  if (left.get_inc() != right.get_inc())
  {
    py_error("Error: left and right part have different incidence media.");
    exit (-1);
  }
  
  uniform = false;

  symmetric = (&left_ex == &right_ex);

  // Determine core.

  materials = left.get_materials();
  vector<Material*> right_mat = right.get_materials();
  materials.insert(materials.end(), right_mat.begin(), right_mat.end());

  unsigned int core_index = 0;
  for (unsigned int i=0; i<materials.size(); i++)
    if ( real(materials[i]->eps_mu()) > real(materials[core_index]->eps_mu()) )
      core_index = i;

  core = materials[core_index];
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::Section2D
//
/////////////////////////////////////////////////////////////////////////////

Section2D::Section2D(const Section2D& section)
{
  params    = section.params;
  left      = section.left;
  right     = section.right;
  core      = section.core;
  uniform   = section.uniform;
  symmetric = section.symmetric;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::operator=
//
/////////////////////////////////////////////////////////////////////////////

Section2D& Section2D::operator=(const Section2D& section)
{
  if (this == &section)
    return *this;

  params    = section.params;
  left      = section.left;
  right     = section.right;
  core      = section.core;
  uniform   = section.uniform;
  symmetric = section.symmetric;
  
  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool Section2D::no_gain_present() const
{
  return left.no_gain_present() && right.no_gain_present();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::eps_at
//
/////////////////////////////////////////////////////////////////////////////

Complex Section2D::eps_at(const Coord& c) const
{
  bool in_left = real(c.c1) < real(left.get_total_thickness());

  Complex c1 = c.c1;
  Limit c1_limit = c.c1_limit;

  if (in_left)
  {
    c1 = left.get_total_thickness() - c.c1;
    c1_limit = (c.c1_limit == Plus) ? Min : Plus;
  }
  else
    c1 -= left.get_total_thickness();

  Coord c_new(c.c2, c.z, c1, c.c2_limit, c.z_limit, c1_limit);

  return in_left ? left.eps_at(c_new) : right.eps_at(c_new);
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::mu_at
//
/////////////////////////////////////////////////////////////////////////////

Complex Section2D::mu_at(const Coord& c) const
{
  bool in_left = c.c1.real() < left.get_total_thickness().real();

  Complex c1 = c.c1;
  Limit c1_limit = c.c1_limit;

  if (in_left)
  {
    c1 = left.get_total_thickness() - c.c1;
    c1_limit = (c.c1_limit == Plus) ? Min : Plus;
  }
  else
    c1 -= left.get_total_thickness();

  Coord c_new(c.c2, c.z, c1, c.c2_limit, c.z_limit, c1_limit);

  return in_left ? left.mu_at(c_new) : right.mu_at(c_new);
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes()
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
// Section2D::find_modes_from_scratch_by_ADR
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes_from_scratch_by_ADR()
{

  cout << "Section2D.find_modes_from_scratch_by_ADR() not yet implemented." 
       << endl;
  
  exit (-1);

#if 0

  // Determine reflection coefficients of walls.

  //SectionWall* l_wall =  leftwall ?  leftwall : global_section. leftwall;
  //SectionWall* r_wall = rightwall ? rightwall : global_section.rightwall;

  //Complex R_left  = l_wall ? l_wall->get_R12() : -1;
  //Complex R_right = r_wall ? r_wall->get_R12() : -1;

  // Find min and max refractive index of structure.

  Complex min_eps_mu = materials[0]->eps()*materials[0]->mu();
  Complex max_eps_mu = materials[0]->eps()*materials[0]->mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu =  materials[i]->eps()*materials[i]->mu();
    
    if (abs(eps_mu) < abs(min_eps_mu))
      min_eps_mu = eps_mu;

    if (abs(eps_mu) > abs(min_eps_mu))
      max_eps_mu = eps_mu;
  }

  // Locate zeros starting from initial contour.

  const Real C = pow(2*pi/global.lambda, 2) / eps0 / mu0;
  
  Real max_beta = abs(sqrt(C*max_eps_mu))+1;

  Real Re=real(global.C_upperright);
  Real Im=imag(global.C_upperright);

  Complex lowerleft (-0.1,      -0.1);
  Complex upperright( Re*max_beta, Im*max_beta);

  SectionDisp disp(left, right, global.lambda, M, symmetric);
  unsigned int zeros = global.N + 2 + materials.size();
  vector<Complex> beta = N_roots(disp, zeros, lowerleft, upperright);
  
  cout << "Calls to slab dispersion relation : "<< disp.times_called() << endl;
  
  // Eliminate false zeros.

  const Real eps_copies = 1e-6;
  
  for (unsigned int i=0; i<n_lossless.size(); i++)
  {
    Complex beta_i = k0*materials[i]->n();
    
    remove_elems(&beta,     beta_i,    eps_copies);
    remove_elems(&beta,    -beta_i,    eps_copies);
    remove_elems(&beta,  Complex(0.0), eps_copies);
  }
  
  // In case of a homogeneous medium, add a TEM mode if appropriate.
  
  //if (materials.size() == 1)
  //  if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
  //                                     && (abs(R_right + 1.0) < 1e-10) )
  //    || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
  //                                     && (abs(R_right - 1.0) < 1e-10) ) )
  //    beta.insert(beta.begin(), 0);
  
  // Check if enough modes are found.

  if (beta.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << beta.size() << "/" << global.N << "). " << endl;
    exit (-1);
  }
  
  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<beta.size(); i++)
  { 
    Complex kz = beta[i];
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Section2D_Mode *newmode 
      = new Section2D_Mode(global.polarisation, kz, this);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();

#endif

}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_from_scratch_by_track
//
/////////////////////////////////////////////////////////////////////////////

// TMP

#include "sectionoverlap.h"

void Section2D::find_modes_from_scratch_by_track()
{
  // Set constants.
  
  const Real eps        = 1e-13;
  const Real eps_zero   = 1e-10;
  const Real eps_copies = 1e-6;
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;
  
  // Determine reflection coefficients of walls.

  //SlabWall* l_wall =  leftwall ?  leftwall : global_slab. leftwall;
  //SlabWall* r_wall = rightwall ? rightwall : global_slab.rightwall;

  //Complex R_left  = l_wall ? l_wall->get_R12() : -1.0;
  //Complex R_right = r_wall ? r_wall->get_R12() : -1.0;

  bool branchcut = false;
  //if ( (abs(R_left) < eps_zero) || (abs(R_right) < eps_zero) )
  //  branchcut = true; // Branchcuts exist only in open structures.

  bool TBC = false; // Transparent boundary conditions.
  //if (   dynamic_cast<SlabWall_TBC*>(l_wall)
  //    || dynamic_cast<SlabWall_TBC*>(r_wall) )
  //  TBC = true;

  bool PC = false; // Photonic crystal boundary conditions.
  //if (   dynamic_cast<SlabWall_PC*>(l_wall)
  //    || dynamic_cast<SlabWall_PC*>(r_wall) )
  //  PC = true;

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

  bool metal = (min_eps_mu_lossless < 0.0);
  
  // Create dispersion relation for lossless structure.
  
  SectionDisp disp(left, right, lambda, M, symmetric);
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

  Complex max_eps_eff;
  if (    (global.polarisation == TM) // Surface plasmon.
       && (max_eps_mu_lossless*min_eps_mu_lossless < 0.0) )
    max_eps_eff = 1.5/(1.0/max_eps_mu_lossless + 1.0/min_eps_mu_lossless);
  else
    max_eps_eff = max_eps_mu_lossless;

  Real prop_kt_end_lossless = abs(sqrt(C*(max_eps_eff - min_eps_mu_lossless)));
  
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
  for (unsigned int i=1; i<=left.get_inc()->N(); i++)
  {
    // TODO: not general for lossy incidence media.

    Complex beta_i = left.get_inc()->get_mode(i)->get_kz();
    Real kt_i = real(C*min_eps_mu_lossless - beta_i*beta_i);

    remove_elems(&kt_prop_lossless,  kt_i, eps_copies);
    remove_elems(&kt_prop_lossless, -kt_i, eps_copies);
  }
  
  // Check if the minima found correspond really to zeros by using a mueller
  // Don't do this when there are branchcuts, as this can cause problems for
  // numerical stability.

  vector<Complex> kt_lossless, kt_complex;
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

      Complex kt_new = mueller(disp, I*kt_prop_lossless[i],
                               I*kt_prop_lossless[i]+0.002,1e-11,0,100,&error);

      if (!error && metal && abs(real(kt_new)) > 0.1)
      {
        py_print("Found complex mode pair.");
        kt_complex.push_back(kt_new);
      }
      else if (!error)
        kt_lossless.push_back(I*imag(kt_new));
    }
  }

  if (kt_lossless.size())
  {
    if (real(sqrt(C*min_eps_mu_lossless - kt_lossless[0]*kt_lossless[0])) > 
      real(left.get_inc()->get_mode(1)->get_kz()) )
    cout << "n_eff higher than expected!" << endl;
    cout << "Slab kz: " << left.get_inc()->get_mode(1)->get_kz() << endl;
  }

  cout << std::setprecision(15);
  cout << "Done lossless guided:" << endl;
  for (unsigned int i=0; i<kt_lossless.size(); i++)
    cout << i << " " << kt_lossless[i] << endl;
  cout << "Calls " << disp.times_called() << endl << std::flush;



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

        if (!error && metal && (abs(imag(kt_new)) > 0.1))
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

  for (unsigned int i=1; i<=left.get_inc()->N(); i++)
  {
    Complex beta_i = left.get_inc()->get_mode(i)->get_kz();
    Complex kt_i = sqrt(C*min_eps_mu_lossless - beta_i*beta_i);

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

  bool degenerate = kt_lossless.size() > kt_lossless_single.size();

  for (unsigned int i=0; i<kt_complex.size(); i++)
  {  
    kt_lossless.push_back(kt_complex[i]);
    kt_lossless.push_back(-real(kt_complex[i])+I*imag(kt_complex[i]));
  }
  
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

  Complex min_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }

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

    for (unsigned int i=1; i<=left.get_inc()->N(); i++)
    {
      Complex beta_i = left.get_inc()->get_mode(i)->get_kz();
      Complex kt_i = sqrt(C*min_eps_mu - beta_i*beta_i);

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
  
  //if (materials.size() == 1)
  //  if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
  //                                     && (abs(R_right + 1.0) < 1e-10) )
  //    || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
  //                                     && (abs(R_right - 1.0) < 1e-10) ) )
  //    kt.insert(kt.begin(), 0);


  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<kt.size(); i++)
  { 
    Complex kz = sqrt(C*min_eps_mu - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    if (real(kt[i]) < -1e-6) // Backward mode.
      kz = -kz;

    disp.set_params(params);
    cout << "n_eff" << kz/2./pi*global.lambda << "kt: " << kt[i] 
         << " f(kt) " << disp(kt[i]) << endl;
    
    Section2D_Mode *newmode = new Section2D_Mode(global.polarisation,kz,this);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }

  sort_modes();
  truncate_N_modes();

  // Test orthogonality.

  return;

  for (unsigned int i=0; i<modeset.size(); i++)
    for (unsigned int j=0; j<modeset.size(); j++)
      if (i !=j)
        std::cout << i << " " << j << " " 
              << overlap(dynamic_cast<Section2D_Mode*>(modeset[i]),
                         dynamic_cast<Section2D_Mode*>(modeset[j]))
                  << std::endl << std::flush;
}





/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_by_sweep
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes_by_sweep()
{   
  // Trace modes from old configuration to new one.
  
  SectionDisp disp(left, right, global.lambda, M, symmetric);
  vector<Complex> params_new = disp.get_params();
    
  vector<Complex> beta_old;
  for (unsigned int i=0; i<modeset.size(); i++)
    beta_old.push_back(modeset[i]->get_kz());
  
  vector<Complex> forbidden;
  vector<Complex> beta =
    traceroot(beta_old,disp,params,params_new,forbidden,global.sweep_steps);

  params = params_new;

  // Check if modes were lost during tracing.

  if (beta.size() < global.N)
  {
    std::ostringstream s;
    s << "Error: didn't find enough modes ("
      << beta.size() << "/" << global.N << ").";
    py_error(s.str());
    exit (-1);
  }

  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<beta.size(); i++)
  {
    Complex kz = beta[i];
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    Section2D_Mode *newmode = new Section2D_Mode(global.polarisation,kz,this);
      
    newmode->normalise();

    modeset.push_back(newmode);
  }
  
  sort_modes();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Section1D::find_modes()
{

  cout << "Section1D.find_modes() not yet implemented." << endl;
  
  exit (-1);

#if 0

  // Determine reflection coefficients of walls.

  SectionWall* l_wall =  leftwall ?  leftwall : global_section. leftwall;
  SectionWall* r_wall = rightwall ? rightwall : global_section.rightwall;

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
  // non-uniform section.

  if (!analytic)
  {
    Section2D tmp = (*core)(get_width()/2.) + (*core)(get_width()/2.);

    if (l_wall)
      tmp.set_left_wall(*l_wall);

    if (r_wall)
      tmp.set_right_wall(*r_wall);
    
    tmp.find_modes();

    for (unsigned int i=1; i<=tmp.N(); i++)
    {
      Section1D_mode *newmode = new 
        Section1D_mode(global.polarisation,tmp.get_mode(i)->get_kz(),this);

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
    Section1D_Mode *newmode = new Section1D_Mode(global.polarisation, k, this);

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
    
    Section1D_Mode *newmode 
      = new Section1D_Mode(global.polarisation, kz, this);

    newmode->normalise();
    
    modeset.push_back(newmode);

    i++;
  }

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;

#endif
}
