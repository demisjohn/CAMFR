
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

#include "section.h"
#include "sectiondisp.h"
#include "sectionmode.h"
#include "sectionoverlap.h"
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

#if 0
  SectionImpl* medium_I  = this;
  SectionImpl* medium_II = dynamic_cast<SectionImpl*>(w);

  // Make sorted list of evaluation points for field cache.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);
  
  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  sort(disc.begin(), disc.end(), RealSorter());

  // Fill field cache.

  const unsigned int N = global.N;

  SectionCache cache(N, disc.size()-1);
  
  for (int i=1; i<=int(N); i++)
  {
    const SectionMode* mode_I
      = dynamic_cast<const SectionMode*>(medium_I ->get_mode(i));
    
    const SectionMode* mode_II
      = dynamic_cast<const SectionMode*>(medium_II->get_mode(i));

    for (int k=0; k<int(disc.size()-1); k++)
    {
      const Coord lower(disc[k],  0,0, Plus);
      const Coord upper(disc[k+1],0,0, Min);

      Complex fw_I_l, bw_I_l, fw_II_l, bw_II_l;
      Complex fw_I_u, bw_I_u, fw_II_u, bw_II_u;

      mode_I ->forw_backw_at(lower, &fw_I_l,  &bw_I_l);    
      mode_I ->forw_backw_at(upper, &fw_I_u,  &bw_I_u);
      mode_II->forw_backw_at(lower, &fw_II_l, &bw_II_l);  
      mode_II->forw_backw_at(upper, &fw_II_u, &bw_II_u);

      cache.fw_l(1,i,k+1) = fw_I_l;
      cache.bw_l(1,i,k+1) = bw_I_l;
      cache.fw_u(1,i,k+1) = fw_I_u;
      cache.bw_u(1,i,k+1) = bw_I_u;

      cache.fw_l(2,i,k+1) = fw_II_l;
      cache.bw_l(2,i,k+1) = bw_II_l;
      cache.fw_u(2,i,k+1) = fw_II_u;
      cache.bw_u(2,i,k+1) = bw_II_u;
    }
  }

  // Calculate overlap matrices

  for (int i=1; i<=int(N); i++)
    for (int j=1; j<=int(N); j++)
    {
      (*O_I_II)(i,j) = overlap
        (dynamic_cast<const SectionMode*>(medium_I ->get_mode(i)),
         dynamic_cast<const SectionMode*>(medium_II->get_mode(j)),
         &cache, &disc, i, j, 1, 2);

      (*O_II_I)(i,j) = overlap
        (dynamic_cast<const SectionMode*>(medium_II->get_mode(i)),
         dynamic_cast<const SectionMode*>(medium_I ->get_mode(j)),
         &cache, &disc, i, j, 2, 1);
      
      if (O_I_I) (*O_I_I)(i,j) = overlap
         (dynamic_cast<const SectionMode*>(medium_I ->get_mode(i)),
          dynamic_cast<const SectionMode*>(medium_I ->get_mode(j)),
          &cache, &disc, i, j, 1, 1);

      if (O_II_II) (*O_II_II)(i,j) = overlap
         (dynamic_cast<const SectionMode*>(medium_II->get_mode(i)),
          dynamic_cast<const SectionMode*>(medium_II->get_mode(j)),
          &cache, &disc, i, j, 2, 2);
    }

#endif
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
      cerr << "Error: expected a slab to initialise a section." << endl;
      exit (-1);
    }    

    uniform = s->is_uniform();
    core = s->get_core();

    return;
  }
  
  // Transparently factor expression into right and left expression.

  cerr << "Automatic splitting of expressions not yet implemented." << endl;
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
    cerr << "Error: left and right part have different incidence media."
         << endl;
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
    if ( real(materials[i]->n()) > real(materials[core_index]->n()) )
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

void Section2D::find_modes_from_scratch_by_track()
{
  // Set constants
  
  const Real eps        = 1e-13;
  const Real eps_zero   = 1e-10;
  const Real eps_copies = 1e-6;
  
  const Real lambda     = global.lambda;
  const Real k0         = 2*pi/lambda;
  
  // Determine reflection coefficients of walls.

  //SectionWall* l_wall =  leftwall ?  leftwall : global_section. leftwall;
  //SectionWall* r_wall = rightwall ? rightwall : global_section.rightwall;

  //Complex R_left  = l_wall ? l_wall->get_R12() : -1;
  //Complex R_right = r_wall ? r_wall->get_R12() : -1;

  bool branchcut = false;
  //if ( (abs(R_left) < eps_zero) || (abs(R_right) < eps_zero) )
  //  branchcut = true; // Branchcuts exist only in open structures.

  bool TBC = false; // Transparent boundary conditions.
  //if (   dynamic_cast<SectionWall_TBC*>(l_wall)
  //    || dynamic_cast<SectionWall_TBC*>(r_wall) )
  //  TBC = true;

  bool PC = false; // Photonic crystal boundary conditions.
  //if (   dynamic_cast<SectionWall_PC*>(l_wall)
  //    || dynamic_cast<SectionWall_PC*>(r_wall) )
  //  PC = true;

  // Find and max/min real refractive indices of lossless structure.

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
  
  SectionDisp disp(left, right, lambda, M, symmetric);
  params = disp.get_params();

  vector<Complex> params_lossless;
  for (unsigned int i=0; i<params.size(); i++)
    params_lossless.push_back(real(params[i]));

  disp.set_params(params_lossless);

  // Make a rough estimate of separation between modes.

  Real d_beta = real(pi/get_width());

  // Find the modes of the structure, proceeding in a number
  // of possible steps:
  //
  //   I:   find guided modes along real axis for lossless structure.
  //   II:  find prop. rad.  modes along real axis for lossless structure.
  //   III: find imag. rad.  along imag axis for lossless structrue.
  //   IV:  track those to modes of true lossy structure.
  //   V:   find modes in complex plane for true structure.


  
  //
  // I: Find guided modes of lossless structure.
  //

  //global.eigen_calc = arnoldi;
  
  Wrap_real_to_abs prop_wrap(disp);

  if (d_beta > k0*max_n)
    d_beta = k0*max_n;

  int sec = (global.precision_enhancement == 1) ? 1 : 0; // security level.

  vector<Real> beta_guided_lossless = brent_all_minima
    (prop_wrap,k0*min_n,k0*max_n,d_beta/global.precision,eps,sec);

  // If precision_enhancement is larger than 1, refine our search for
  // the guided modes.

  if (global.precision_enhancement > 1)
  {
    vector<Real> beta_fine = brent_refine_minima
      (prop_wrap,beta_guided_lossless,global.dx_enhanced,
       d_beta/global.precision/global.precision_enhancement,eps,1);
    beta_guided_lossless = beta_fine;
  }
  
  reverse(beta_guided_lossless.begin(), beta_guided_lossless.end());

  // Eliminate false zeros.

  remove_elems(&beta_guided_lossless, 0.0, eps_copies);

  global_slab.beta = 0.0;
  for (unsigned int i=1; i<=left.get_inc()->N(); i++)
  {    
    Real beta_i = real(left.get_inc()->get_mode(i)->get_kz());

    remove_elems(&beta_guided_lossless,  beta_i, eps_copies);
    remove_elems(&beta_guided_lossless, -beta_i, eps_copies);
  }

  if (beta_guided_lossless[0] > 
      real(left.get_inc()->get_mode(1)->get_kz()))
    cout << "n_eff higher than expected!" << endl;
  cout << "Slab kz: " << left.get_inc()->get_mode(1)->get_kz() << endl;

  cout << setprecision(15);
  cout << "Done lossless:" << endl;
  for (unsigned int i=0; i<beta_guided_lossless.size(); i++)
    cout << i << beta_guided_lossless[i]/2/pi*global.lambda << endl;
  cout << "Calls " << disp.times_called() << endl;
  
  // Check if the minima found correspond really to zeros and increase
  // precision using a mueller solver (in the absense of branchcuts).
  // Don't do this if precision_enhancement is larger than 1, since this
  // typically indicate the presence of nearly degenerate modes.

  vector<Complex> beta_lossless;
  for (unsigned int i=0; i<beta_guided_lossless.size(); i++)
  { 
    const Real fx = abs(disp(beta_guided_lossless[i]));
    
    if (fx > 1e-2)
    {
      cout << "Warning: possibly insufficient precision around beta "
           << beta_guided_lossless[i] << "." << endl;
      cout << "Removing this candidate with abs(fx) " << fx << "." << endl;
    }
    else
    {
      Real beta_new;
      
      if ( (!branchcut) && (global.precision_enhancement == 1) )
        beta_new = real(mueller(disp, beta_guided_lossless[i],
                                beta_guided_lossless[i]+0.002*I));
      else
        beta_new = beta_guided_lossless[i];

      beta_lossless.push_back(beta_new);
    }
  }

  cout << "Done refining" << endl;
  for (unsigned int i=0; i<beta_lossless.size(); i++)
    cout << i << " " << beta_lossless[i] /2. /pi * global.lambda << endl;
  cout << "Calls " << disp.times_called() << endl;

  //
  // II: Find propagating radiation modes of lossless structure.
  //

  global.eigen_calc = lapack;

  vector<Real> beta_prop_rad_lossless;
  if (beta_lossless.size() < global.N)
    beta_prop_rad_lossless = brent_all_minima
     (prop_wrap,0.0001,k0*min_n,d_beta/global.precision_rad,eps,sec);

  reverse(beta_prop_rad_lossless.begin(), beta_prop_rad_lossless.end());

  cout << "Done prop rad" << endl;
      for (unsigned int i=0; i<beta_prop_rad_lossless.size(); i++)
        cout << beta_prop_rad_lossless[i] /2./pi*global.lambda << endl;
  cout << "Calls prop rad "  << disp.times_called() << endl;
  
  // Check if the minima found correspond really to zeros and increase
  // precision using a mueller solver (in the absense of branchcuts).

  for (unsigned int i=0; i<beta_prop_rad_lossless.size(); i++)
  { 
    const Real fx = abs(disp(beta_prop_rad_lossless[i]));
    
    if (fx > 1e-2)
    {
      cout << "Warning: possibly insufficient precision around beta "
           << beta_prop_rad_lossless[i] << "." << endl;
      cout << "Removing this candidate with abs(fx) " << fx << "." << endl;
    }
    else
    {
      Real beta_new;
      
      if (!branchcut)
        beta_new = real(mueller(disp, beta_prop_rad_lossless[i],
                                beta_prop_rad_lossless[i]+0.002*I));
      else
        beta_new = beta_prop_rad_lossless[i];

      beta_lossless.push_back(beta_new);
    }
  }


  cout << "Done prop rad refin" << endl;
  for (unsigned int i=0; i<beta_lossless.size(); i++)
    cout << beta_lossless[i] /2./pi*global.lambda << endl;
  cout << "Calls prop rad refined " << disp.times_called() << endl;

  
  //
  // III: Find evanescent radiation modes of lossless structure.
  //

  if (!(branchcut || TBC))
  {
    Wrap_imag_to_abs evan_wrap(disp);
    
    Real beta_begin = .0001;
    int extra = 1;
    int modes_left = (beta_lossless.size() >= global.N) ? 0
      : global.N - beta_lossless.size();
    
    while (modes_left)
    {
      if (extra > 1)
        cout << "Lost too many candidate modes. Expanding search region."
             << endl;
      
      vector<Real> beta_evan_lossless = brent_N_minima
        (evan_wrap,beta_begin,modes_left+extra,
         d_beta/global.precision_rad,eps,1);

      cout << "Done evan " << endl;
      for (unsigned int i=0; i<beta_evan_lossless.size(); i++)
        cout << beta_evan_lossless[i] /2./pi*global.lambda << endl;
      cout << "Calls evan rad " << disp.times_called() << endl;
      
      // Check if the minima found correspond really to zeros and increase
      // precision with mueller solver.

      for (unsigned int i=0; i<beta_evan_lossless.size(); i++)
      { 
        const Real fx = abs(disp(I*beta_evan_lossless[i]));
        if (fx > 1e-2)
        {
          cout << "Warning: possibly insufficient precision around beta "
               <<  beta_evan_lossless[i] << "." << endl;
          cout << "Removing this candidate with abs(fx) "
               << fx << "." << endl;
        }
        else
        {
          Real beta_new = branchcut
            ? beta_evan_lossless[i]
            : imag(mueller(disp,I*beta_evan_lossless[i],
                           I*beta_evan_lossless[i]+0.002));

          beta_lossless.push_back(I*beta_new);
        }
      }

      modes_left = (beta_lossless.size() >= global.N + 2) ? 0
        : global.N - beta_lossless.size() + 2;

      beta_begin = real(beta_lossless[beta_lossless.size()-1])+.0001;
      extra++;
    }
  }
 
  cout << "Done evan refin" << endl;
  for (unsigned int i=0; i<beta_lossless.size(); i++)
    cout << i << " " << beta_lossless[i] /2./pi*global.lambda 
         << " " << beta_lossless[i] << endl;
  cout << "Calls evan refined " << disp.times_called() << endl;

  
  
  //
  // IV: Eliminate doubles and false zeros and trace modes to those of
  //      true structure.
  //

  if (global.precision_enhancement == 1)
    remove_copies(&beta_lossless, eps_copies);

  remove_elems(&beta_lossless, Complex(0.0), eps_copies);

  if (beta_lossless.size() > global.N)
    beta_lossless.erase(beta_lossless.begin()+global.N, beta_lossless.end());

  global_slab.beta = 0.0;
  vector<Complex> forbidden;
  for (unsigned int i=1; i<=left.get_inc()->N(); i++)
  {
    Complex beta_i = left.get_inc()->get_mode(i)->get_kz();

    remove_elems(&beta_lossless,  beta_i, eps_copies);
    remove_elems(&beta_lossless, -beta_i, eps_copies);

    forbidden.push_back( beta_i);
    forbidden.push_back(-beta_i);
  }
  
  vector<Complex> beta;
  if (global.chunk_tracing == true)
    beta = traceroot_chunks
      (beta_lossless,disp,params_lossless,params,forbidden,global.sweep_steps);
  else
    beta = traceroot
      (beta_lossless,disp,params_lossless,params,forbidden,global.sweep_steps);

  cout << "Calls traceroot " << disp.times_called() << endl;
 

  //
  // V: Find modes in complex plane
  //

#if 0

  disp.set_params(params);
  
  if ( (branchcut || TBC) && (beta.size() < global.N) )
  {
    Real evan_beta_end = 2*(global.N-beta.size())*d_beta; // Not general!

    Real Re=real(global.C_upperright);
    Real Im=imag(global.C_upperright);

    Complex lowerleft (0.001,          0.001);
    //Complex upperright(Re*evan_beta_end, Im*k0*max_n);
    Complex upperright = global.C_upperright;

    const int sections = int(global.C_steps);

    vector<Complex> beta_complex = allroots(disp, lowerleft, upperright);

    // Eliminate doubles and false zeros.

    remove_copies(&beta_complex, eps_copies);

    for (unsigned int i=0; i<n_lossless.size(); i++)
    {
      Complex beta_i = k0*n_lossless[i]);

      remove_elems(&beta_complex,  beta_i, eps_copies);
      remove_elems(&beta_complex, -beta_i, eps_copies);
    }

    cout << "Found " << beta_complex.size() << " complex zeros in region "
         << lowerleft << " " << upperright << "." <<endl;

    for (unsigned int i=0; i<beta_complex.size(); i++)
      beta.push_back(beta_complex[i]);
  }
#endif



  //
  // Finalise.
  //
  
  // In case of a homogeneous medium, add a TEM mode if appropriate.
/*  
  if (materials.size() == 1)
    if ( ( (global.polarisation == TM) && (abs(R_left  + 1.0) < 1e-10)
                                       && (abs(R_right + 1.0) < 1e-10) )
      || ( (global.polarisation == TE) && (abs(R_left  - 1.0) < 1e-10)
                                       && (abs(R_right - 1.0) < 1e-10) ) )
      beta.insert(beta.begin(), 0);
*/  
  // Check if enough modes are found.

  if (beta.size() < global.N)
  {
    cout << "Error: didn't find enough modes ("
         << beta.size() << "/" << global.N << "). " << endl;
    //exit (-1);
  }
  
  // Create modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();
  
  for (unsigned int i=0; i<beta.size(); i++)
  { 
    Complex kz = beta[i];

    disp.set_params(params);
    cout << "kz: " << kz << " f(kz) " << disp(kz) << endl;
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;
    
    cVector fw_field(global.N, fortranArray);
    
    Section2D_Mode *newmode 
      = new Section2D_Mode(global.polarisation, kz, this, fw_field);
    
    newmode->normalise();
    
    modeset.push_back(newmode);
  }
  
  sort_modes();
  truncate_N_modes();
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

    cVector fw_field(global.N, fortranArray);
    
    Section2D_Mode *newmode 
      = new Section2D_Mode(global.polarisation, kz, this, fw_field);
      
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
