
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
#include <algorithm>
#include "section.h"
#include "refsection.h"
#include "sectiondisp.h"
#include "sectionmode.h"
#include "sectionoverlap.h"
#include "../slab/slabmatrixcache.h"
#include "../slab/isoslab/slaboverlap.h"
#include "../slab/isoslab/slabmode.h"

using std::vector;
using std::cout;
using std::endl;

#include "../../util/vectorutil.h"
#include "../../util/index.h"

/////////////////////////////////////////////////////////////////////////////
//
// SectionGlobal
//
/////////////////////////////////////////////////////////////////////////////

SectionGlobal global_section = {0.0, 0.0, E_wall, E_wall, false};



/////////////////////////////////////////////////////////////////////////////
//
// SectionImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

void SectionImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I, cMatrix* O_II_II)
{
  Section2D* medium_I  = dynamic_cast<Section2D*>(this);
  Section2D* medium_II = dynamic_cast<Section2D*>(w);
  
  // Widths equal?

  const Real eps = 1e-10; // Don't choose too low.
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << medium_I ->get_width() << " and " << medium_II->get_width();
    py_error(s.str());
    return;
  }

  // Make sorted list of separation points between different slabs.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);

  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  sort(disc.begin(), disc.end(), RealSorter());

  // Calculate overlap matrices.

  *O_I_II = 0.0;  
  *O_II_I = 0.0;

  if (O_I_I) 
    *O_I_I = 0.0;
  if (O_II_II) 
    *O_II_II = 0.0;
  
  for (unsigned int k=0; k<disc.size()-1; k++) // Loop over slices.
  {
    // Create field cache.

    vector<FieldExpansion> field_I, field_II;
    
    for (int i=1; i<=int(global.N); i++)
    {
      cVector fw_I(medium_I->get_M2(),fortranArray);
      cVector bw_I(medium_I->get_M2(),fortranArray);

      dynamic_cast<SectionMode*>(medium_I->get_mode(i))
        ->get_fw_bw(disc[k], Plus, &fw_I, &bw_I);

      field_I.push_back(FieldExpansion(NULL, fw_I, bw_I));

      cVector fw_II(medium_II->get_M2(),fortranArray);
      cVector bw_II(medium_II->get_M2(),fortranArray);

      dynamic_cast<SectionMode*>(medium_II->get_mode(i))
        ->get_fw_bw(disc[k], Plus, &fw_II, &bw_II);

      field_II.push_back(FieldExpansion(NULL, fw_II, bw_II));
    }

    // Create overlap cache.

    if (medium_I->get_M2() != medium_II->get_M2())
    {
      std::cerr 
        << "Error: cache not yet general enough to deal with different M."
        << std::endl;
      return;
    }

    SlabImpl* slab_I = medium_I->
      slabs[index_lookup(disc[k],Plus,medium_I ->discontinuities)]->get_impl();

    SlabImpl* slab_II = medium_II->
      slabs[index_lookup(disc[k],Plus,medium_II->discontinuities)]->get_impl();
    
    vector<Complex> disc_slab = slab_I->disc_intersect(slab_II);
    
    SlabCache cache(medium_I->get_M2(), disc_slab.size()-1);
    slab_I->fill_field_cache(&cache, slab_II, disc_slab);

    OverlapMatrices m(slab_I, slab_II, &cache, &disc_slab, false);

    // Calc overlap slice.

    for (int i=1; i<=int(global.N); i++)   
      for (int j=1; j<=int(global.N); j++)
      {
        (*O_I_II)(i,j) += overlap_slice
          (dynamic_cast<SectionMode*>(medium_I ->get_mode(i)),
           dynamic_cast<SectionMode*>(medium_II->get_mode(j)),
           disc[k], disc[k+1], &field_I[i-1], &field_II[j-1], &m, 1, 2);

        (*O_II_I)(i,j) += overlap_slice
          (dynamic_cast<SectionMode*>(medium_II->get_mode(i)),
           dynamic_cast<SectionMode*>(medium_I ->get_mode(j)),
           disc[k], disc[k+1], &field_II[i-1], &field_I[j-1], &m, 2, 1);

        if (O_I_I) (*O_I_I)(i,j) += overlap_slice
          (dynamic_cast<SectionMode*>(medium_I ->get_mode(i)),
           dynamic_cast<SectionMode*>(medium_I ->get_mode(j)),
           disc[k], disc[k+1], &field_I[i-1], &field_I[j-1]);
      
        if (O_II_II) (*O_II_II)(i,j) += overlap_slice
         (dynamic_cast<SectionMode*>(medium_II->get_mode(i)),
          dynamic_cast<SectionMode*>(medium_II->get_mode(j)),
          disc[k], disc[k+1], &field_II[i-1], &field_II[j-1]);
      }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(const Term& t) : leftwall_sc(NULL), rightwall_sc(NULL)
{
  Slab* slab = dynamic_cast<Slab*>(t.get_wg());
  Complex d = t.get_d()
    + I*global_section.left_PML + I*global_section.right_PML;

  if (slab)
    s = new Section1D(*slab, d);
  else
  {
    py_error("Error: expected a slab to initialise a section.");
    return;
  }

  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(Expression& expression, int M1, int M2)
  : leftwall_sc(NULL), rightwall_sc(NULL)
{
  //
  // 1D section
  //

  if (expression.get_size() == 1)
  {
    Slab* slab = dynamic_cast<Slab*>(expression.get_term(0)->get_wg());
    Complex d = expression.get_term(0)->get_d() 
      + I*global_section.left_PML + I*global_section.right_PML;

    if (slab)
      s = new Section1D(*slab, d);
    else
    {
      py_error("Error: expected a slab to initialise a section.");
      return;
    }    

    uniform = s->is_uniform();
    core = s->get_core();

    return;
  }

  //
  // 2D section.
  //
  
  // Find core.

  Expression ex = expression.flatten();

  Complex max_eps = 0.0;
  unsigned int max_eps_i = 0;
  for (unsigned int i=1; i<ex.get_size(); i++)
  {
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());
    if (s)
    {
      Complex eps_i = s->eps_avg();

      if (real(eps_i) > real(max_eps))
      {
        max_eps = eps_i;
        max_eps_i = i;
      }
    }
  }

  if (max_eps_i == 1)
    max_eps_i += 2;

  // Create right hand side expression.

  Expression right_ex;
  for (int i=max_eps_i; i<ex.get_size(); i++)
  {
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());

    if (s)
    {
      Complex d = ex.get_term(i)->get_d();

      if (i == ex.get_size()-1)
        d += I*global_section.right_PML;

      right_ex += Term((*s)(d));
    }
  }

  if (global_section.rightwall == E_wall)
    rightwall_sc = new E_Wall(*right_ex.get_ext());
  if (global_section.rightwall == H_wall)
    rightwall_sc = new H_Wall(*right_ex.get_ext());
  right_ex.add_term(Term(*rightwall_sc));

  // Create left hand side expression.

  Expression left_ex;

  Slab* slab = dynamic_cast<Slab*>(ex.get_term(max_eps_i)->get_wg());
  Scatterer* sc 
    = interface_cache.get_interface(slab,slab);
  left_ex.add_term(Term(*sc));

  left_ex += Term((*slab)(0.0));

  for (int i=max_eps_i-1; i>=0; i--)
  { 
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());

    if (s)
    {
      Complex d = ex.get_term(i)->get_d();

      if (i == 1)
        d += I*global_section.left_PML;

      left_ex += Term((*s)(d));
    }
  }

  left_ex.remove_term_front();
  left_ex.remove_term_front();

  if (global_section.leftwall == E_wall)
    leftwall_sc = new E_Wall(*left_ex.get_ext());
  if (global_section.leftwall == H_wall)
    leftwall_sc = new H_Wall(*left_ex.get_ext());
  left_ex.add_term(Term(*leftwall_sc));  

  // Create Section.

  s = new Section2D(left_ex, right_ex, M1, M2);

  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(Expression& left_ex, Expression& right_ex, int M1, int M2)
{
  bool symmetric = (&left_ex == &right_ex);

  // Add walls.

  if (global_section.leftwall == E_wall)
    leftwall_sc = new E_Wall(*left_ex.get_ext());
  if (global_section.leftwall == H_wall)
    leftwall_sc = new H_Wall(*left_ex.get_ext());  
  left_ex.add_term(Term(*leftwall_sc));

  if (!symmetric)
  {
    if (global_section.rightwall == E_wall)
      rightwall_sc = new E_Wall(*right_ex.get_ext());
    if (global_section.rightwall == H_wall)
      rightwall_sc = new H_Wall(*right_ex.get_ext());  
    right_ex.add_term(Term(*rightwall_sc));
  }
  else
    rightwall_sc = 0;

  // Create Section2D

  if (    (abs(global_section. left_PML) > 1e-6) 
       || (abs(global_section.right_PML) > 1e-6) ) 
    std::cout 
      << "Warning: left or right PML not supported from this constructor." 
      << std::endl;

  s = symmetric ? new Section2D(left_ex, left_ex,  M1, M2) 
                : new Section2D(left_ex, right_ex, M1, M2);

  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::Section2D
//
/////////////////////////////////////////////////////////////////////////////

Section2D::Section2D
  (Expression& left_ex, Expression& right_ex, int M1_, int M2_)
    : left(left_ex), right(right_ex)
{
  // Check values.

  if (left.get_inc() != right.get_inc())
  {
    py_error("Error: left and right part have different incidence media.");
    return;
  }

  M1 = M1_;
  M2 = M2_;

  uniform = false;

  symmetric =     (&left_ex == &right_ex) 
              && (global_section.leftwall == global_section.rightwall);

  // Determine core.

  materials = left.get_materials();
  vector<Material*> right_mat = right.get_materials();
  materials.insert(materials.end(), right_mat.begin(), right_mat.end());

  unsigned int core_index = 0;
  for (unsigned int i=0; i<materials.size(); i++)
    if ( real(materials[i]->eps_mu()) > real(materials[core_index]->eps_mu()) )
      core_index = i;

  core = materials[core_index];

  // Determine slabs.

  Complex z = 0; //I*global_section.left_PML;

  Expression left_flat = left_ex.flatten();
  for (int i=left_flat.get_size()-1; i>=0; i--)
  {
    const Term* t = left_flat.get_term(i);
    if (t->get_type() != WAVEGUIDE)
      continue;
    
    z += t->get_d();

    discontinuities.push_back(z);
    slabs.push_back(dynamic_cast<Slab*>(t->get_wg()));
  }
  
  Expression right_flat = right_ex.flatten();
  for (unsigned int i=0; i<right_flat.get_size(); i++)
  {
    const Term* t = right_flat.get_term(i);
    if (t->get_type() != WAVEGUIDE)
      continue;    

    z += t->get_d();

    discontinuities.push_back(z);
    slabs.push_back(dynamic_cast<Slab*>(t->get_wg()));
  }

  Complex z_back = discontinuities.back(); // + I*global_section.right_PML;
  discontinuities.pop_back();
  discontinuities.push_back(z_back);
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
  {
    if (global.solver == series)
      find_modes_from_series();
    else
      find_modes_from_scratch_by_track();
  }

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_from_series
//
/////////////////////////////////////////////////////////////////////////////
  
struct sorter
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
      {return ( real(beta_a) > real(beta_b) );}
};

void Section2D::find_modes_from_series()
{
  // Check M1.

  int n = M1/2;
  
  if (2*n != M1)
    py_print("Warning: changing M1 to even number.");

  M1 = 2*n;

  // Find min eps mu.

  Complex min_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }

  const Real C0 = pow(2*pi/global.lambda, 2) / eps0 / mu0;

  // Calc overlap matrices.

  Material ref_mat(1.0);
  RefSection ref(ref_mat, get_width(), get_height(), M1);

  cMatrix O_EE(n,n,fortranArray); cMatrix O_MM(n,n,fortranArray);
  cMatrix O_EM(n,n,fortranArray); cMatrix O_zz(n,n,fortranArray);

  ref.calc_overlap_matrices(this, &O_EE, &O_MM, &O_EM, &O_zz);

  // Form auxiliary matrix.

  const Real omega = 2*pi/global.lambda * c;

  cMatrix inv_O_zz(n,n,fortranArray);
  if(global.stability != SVD)
    inv_O_zz.reference(invert    (O_zz));
  else
    inv_O_zz.reference(invert_svd(O_zz));

  cMatrix T(n,n,fortranArray);
  for (int i=1; i<=n; i++)
  {
    RefSectionMode* TM_i = dynamic_cast<RefSectionMode*>(ref.get_mode(n+i));

    for (int j=1; j<=n; j++)
    {
      RefSectionMode* TM_j = dynamic_cast<RefSectionMode*>(ref.get_mode(n+j));
      T(i,j) =   TM_i->kt2()/TM_i->get_kz() * inv_O_zz(i,j) 
               * TM_j->kt2()/TM_j->get_kz();
    }

    T(i,i) += pow(omega,3) * ref_mat.eps() * ref_mat.mu() / TM_i->get_kz();
  }

  // Create submatrices for eigenvalue problem.

  cMatrix A(n,n,fortranArray); cMatrix B(n,n,fortranArray);
  cMatrix C(n,n,fortranArray); cMatrix D(n,n,fortranArray);

  A.reference(multiply(T, O_MM));
  
  cMatrix trans_O_EM(n,n,fortranArray);
  trans_O_EM.reference(transpose(O_EM));
  B.reference(multiply(T, trans_O_EM));

  C = O_EM;
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      C(i,j) *= omega * ref.get_mode(i)->get_kz();
  
  D = O_EE;
  for (int i=1; i<=n; i++)
  {
    RefSectionMode* TE_i = dynamic_cast<RefSectionMode*>(ref.get_mode(i));
   
    for (int j=1; j<=n; j++)
      D(i,j) *= omega * TE_i->get_kz();

    D(i,i) -= TE_i->kt2();
  }

  // Solve eigenvalue problem.

  cMatrix E(M1,M1,fortranArray);
  
  blitz::Range r1(1,n); blitz::Range r2(n+1,2*n);

  E(r1,r1) = A; E(r1,r2) = B;
  E(r2,r1) = C; E(r2,r2) = D;

  cVector kz2(M1,fortranArray);
  kz2 = eigenvalues(E);

  // Refine estimates.

  vector<Complex> kz2_coarse;
  for (unsigned int i=1; i<=M1; i++)
    kz2_coarse.push_back(kz2(i));

  std::sort(kz2_coarse.begin(), kz2_coarse.end(),sorter());

  for (unsigned int i=0; i<kz2_coarse.size(); i++)
    std::cout << "coarse" << i << " " << kz2_coarse[i] << std::endl;

  kz2_coarse.erase(kz2_coarse.begin()+global.N, kz2_coarse.end());

  vector<Complex> kt_coarse;
  for (unsigned int i=0; i<kz2_coarse.size(); i++)
  {
    Complex kt = sqrt(C0*min_eps_mu - kz2_coarse[i]);

    if (imag(kt) < 0) 
      kt = -kt;

    if (abs(imag(kt)) < 1e-12)
      if (real(kt) > 0)
        kt = -kt;

    kt_coarse.push_back(kt);
  }

  SectionDisp disp(left, right, global.lambda, M2, symmetric);
  vector<Complex> kt = mueller(disp, kt_coarse, 1e-8, 50);

  // Eliminate false zeros.

  for (unsigned int i=1; i<=left.get_inc()->N(); i++)
  {
    Complex beta_i = left.get_inc()->get_mode(i)->get_kz();
    Complex kt_i = sqrt(C0*min_eps_mu - beta_i*beta_i);

    remove_elems(&kt,  kt_i, 1e-6);
    remove_elems(&kt, -kt_i, 1e-6);
  }

  // Create modeset.

  for (unsigned int i=0; i<kt.size(); i++)
  { 
    Complex kz = sqrt(C0*min_eps_mu - kt[i]*kt[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    Section2D_Mode* newmode
     = new Section2D_Mode(global.polarisation,kz,this);

    newmode->normalise();

    modeset.push_back(newmode);
  }

  sort_modes();
  truncate_N_modes();  
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
  
  SectionDisp disp(left, right, lambda, M2, symmetric);
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

      if (abs(disp(kt_new)) > 1e-2)
        error = true;

      if (!error && abs(real(kt_new)) > 0.001)
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
      real(left.get_inc()->get_mode(1)->get_kz())+0.1 )
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

  if (!(branchcut || TBC || global_section.guided_only))
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

        std::cout << kt_new << " -> " << abs(disp(kt_new)) << std::endl;
        
        if (abs(disp(kt_new)) > 1e-2)
        {
          std::cout << "rejected" << std::endl;
          error = true;
        }

        if (!error && (abs(imag(kt_new)) > 0.001))
        {
          py_print("Found complex mode pair.");
          kt_complex.push_back(kt_new);
        }
        else if (!error)
          kt_lossless.push_back(real(kt_new));
      }

      int modes_needed = global.N + 2*left.get_inc()->N() + 2;
      
      modes_left = (kt_lossless.size() + 2*kt_complex.size() >= modes_needed) 
        ? 0
        : modes_needed - kt_lossless.size() -2*kt_complex.size();

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

  disp.set_params(params);
  
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

    cout << "n_eff" << kz/2./pi*global.lambda << "kt: " << kt[i] 
         << " f(kt) " << disp(kt[i]) << endl;
    
    Section2D_Mode *newmode = new Section2D_Mode(global.polarisation,kz,this);
    
    newmode->normalise();   
    modeset.push_back(newmode);
  }

  sort_modes();
  truncate_N_modes();

  // TODO: check if N modes are found.

  // Test orthogonality.
/*
  return;

  cMatrix O(modeset.size(), modeset.size(), fortranArray);
  overlap_matrices(&O, this, this);
  std::cout << O << std::endl;

  return;

  for (unsigned int i=0; i<modeset.size(); i++)
    for (unsigned int j=0; j<modeset.size(); j++)
        std::cout << i << " " << j << " " 
              << overlap(dynamic_cast<Section2D_Mode*>(modeset[i]),
                         dynamic_cast<Section2D_Mode*>(modeset[j]))
                  << std::endl << std::flush;

  exit(-1);
*/
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_by_sweep
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes_by_sweep()
{   
  // Trace modes from old configuration to new one.
  
  SectionDisp disp(left, right, global.lambda, M2, symmetric);
  vector<Complex> params_new = disp.get_params();
    
  vector<Complex> beta_old;
  for (unsigned int i=0; i<modeset.size(); i++)
    beta_old.push_back(modeset[i]->get_kz());
  
  vector<Complex> forbidden;
  vector<Complex> beta =
    traceroot(beta_old,disp,params,params_new,forbidden,global.sweep_steps);

  params = params_new;

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

  truncate_N_modes();

  // Gracefully handle case when not enough modes were found.

  if (modeset.size() < global.N)
  {
    std::ostringstream s;
    s << "Error: didn't find enough modes ("
      << modeset.size() << "/" << global.N << "). ";
    py_error(s.str());
    
    while (modeset.size() != global.N)
    {
      Section2D_Mode *dummy = new Section2D_Mode(unknown, 0.0, this);
      modeset.push_back(dummy);
    }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Section1D::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Section1D::find_modes()
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

  slab->find_modes();

  // Determine reflection coefficients of walls.

  Complex R_left  = global_section. leftwall == E_wall ? -1.0 : 1.0;
  Complex R_right = global_section.rightwall == E_wall ? -1.0 : 1.0;

  // Clear modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Add modes.
  
  const Complex start  = (abs(R_left*R_right-1.0) < 1e-10) ? 0 : pi;

  for (unsigned int i=1; i<=global.N; i++)
  {
    SlabMode* mode = dynamic_cast<SlabMode*>(slab->get_mode(i));

    for (unsigned int j=0; j<=global.N; j++)
    {
      Complex kz = sqrt(pow(mode->get_kz0(),2) - pow((start+2*j*pi)/2./d,2));

      if (real(kz) < 0) 
        kz = -kz;

      if (abs(real(kz)) < 1e-12)
        if (imag(kz) > 0)
          kz = -kz;
    
      if (abs(kz-2*pi/global.lambda*core->n()) > 1e-6)
      {
        Section1D_Mode *newmode
          = new Section1D_Mode(slab->get_mode(i)->pol, kz, mode, this);

        newmode->normalise();
    
        modeset.push_back(newmode);
      }
    }
  }

  sort_modes();
  truncate_N_modes();

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}
