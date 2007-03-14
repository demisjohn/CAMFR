
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
#include <fstream>
#include <algorithm>
#include "section.h"
#include "refsection.h"
#include "sectiondisp.h"
#include "sectionmode.h"
#include "sectionoverlap.h"
#include "../slab/slabmatrixcache.h"
#include "../slab/isoslab/slabmode.h"
#include "../slab/isoslab/slaboverlap.h"
#include "../../math/calculus/fourier/fourier.h"

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

SectionGlobal global_section = {0.0,0.0,E_wall,E_wall,L,none,0,0,2.0,
   false,0.5,1.0,false,true,false,true,true,false,10.,50,0.0,0.0,1.0,false};


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

  // Dimensions equal?

  const Real eps = 1e-10; // Don't choose too low.
  
  if (abs(medium_I->get_width() - medium_II->get_width()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << medium_I ->get_width() << " and " << medium_II->get_width();
    py_error(s.str());
    return;
  }

  if (abs(medium_I->get_height() - medium_II->get_height()) > eps)
  {
    std::ostringstream s;
    s << "Warning: complex heights don't match: "
      << medium_I ->get_height() << " and " << medium_II->get_height();
    py_error(s.str());
    return;
  }

  if (abs(medium_I->N() - medium_II->N()) > eps)
  {
    std::ostringstream s;
    s << "Error: diffent number of modes in sections: "
      << medium_I->N() << " and " << medium_II->N();
    py_error(s.str());
    return;
  }
 
  // Overlap for partial mode correction case.

  if (global_section.mode_correction != full)
  {
    *O_I_II = 0.0;  
    *O_II_I = 0.0;    

    if (O_I_I) 
      *O_I_I = 0.0;
    if (O_II_II) 
      *O_II_II = 0.0;

    for (int i=1; i<=int(medium_I->N()); i++)   
      for (int j=1; j<=int(medium_I->N()); j++)
      { 
        (*O_I_II)(i,j) = overlap
          (dynamic_cast<Section2D_Mode*>(medium_I ->get_mode(i)),
           dynamic_cast<Section2D_Mode*>(medium_II->get_mode(j)));

        (*O_II_I)(i,j) = overlap
          (dynamic_cast<Section2D_Mode*>(medium_II->get_mode(i)),
           dynamic_cast<Section2D_Mode*>(medium_I ->get_mode(j)));
        
        if (O_I_I) (*O_I_I)(i,j) = overlap
          (dynamic_cast<Section2D_Mode*>(medium_I ->get_mode(i)),
           dynamic_cast<Section2D_Mode*>(medium_I ->get_mode(j)));
      
        if (O_II_II) (*O_II_II)(i,j) = overlap
         (dynamic_cast<Section2D_Mode*>(medium_II->get_mode(i)),
          dynamic_cast<Section2D_Mode*>(medium_II->get_mode(j)));
      }

    py_print("Done overlap.");
    return;
  }

  // Make sorted list of separation points between different slabs.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);

  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  std::sort(disc.begin(), disc.end(), RealSorter());

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
  //if ( (M1 < M2) || (M2 < global.N) )
  //  py_print("Warning: M1 > M2 > N not fulfilled.");

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

  global.always_dense = true; // Uniform slabs are not orthogonal here.
  
  // Find core.

  Expression ex = expression.flatten();

  Complex max_eps = -1e9;
  unsigned int max_eps_i = 0;

  int i_start = 1;
  int i_stop = ex.get_size()-1;

  // Make sure there are no uniform slabs at the edges that end up as
  // core, otherwise they will have a zero reflection matrix and the
  // mode finding will fail.

  for (unsigned int i=i_start; i<=i_stop; i++)
  {
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());
    if (s && s->is_uniform() && global_section.leftwall == no_wall)
    {
      i_start++;
      break;
    }
  }
  
  for (unsigned int i=i_stop; i>=i_start; i--)
  {
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());
    if (s && s->is_uniform() && global_section.rightwall == no_wall)
    {
      i_stop--;
      break;
    }
  }

  // Find actual core.
  
  for (unsigned int i=i_start; i<=i_stop; i++)
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
  
  // Create right hand side expression.

  Expression right_ex;
  for (int i=max_eps_i; i<ex.get_size(); i++)
  {
    Slab* s = dynamic_cast<Slab*>(ex.get_term(i)->get_wg());

    if (s)
    {
      Complex d = ex.get_term(i)->get_d();

      if (i==max_eps_i)
        d /= 2.;

      if (i == ex.get_size()-1)
        d += I*global_section.right_PML;

      right_ex += Term((*s)(d));
    }
  }

  if (global_section.rightwall == E_wall)
    rightwall_sc = new E_Wall(*right_ex.get_ext());
  if (global_section.rightwall == H_wall)
    rightwall_sc = new H_Wall(*right_ex.get_ext());
  if (global_section.rightwall != no_wall)
    right_ex.add_term(Term(*rightwall_sc));

  // Create left hand side expression.

  Expression left_ex;

  Slab* slab0 = dynamic_cast<Slab*>(ex.get_term(max_eps_i)->get_wg());
  Complex d0 = ex.get_term(max_eps_i)->get_d()/2.;

  if (max_eps_i == 1)
    d0 += I*global_section.left_PML;

  left_ex += Term((*slab0)(d0));

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

  if (global_section.leftwall == E_wall)
    leftwall_sc = new E_Wall(*left_ex.get_ext());
  if (global_section.leftwall == H_wall)
    leftwall_sc = new H_Wall(*left_ex.get_ext());
  if (global_section.leftwall != no_wall)
    left_ex.add_term(Term(*leftwall_sc));

  // Create Section.

  s = new Section2D(left_ex, right_ex, M1, M2);
  uniform = s->is_uniform();
  core = s->get_core();

  global.always_dense = false;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section::Section
//
//  There are situations where the core finding algorithm from the
//  previous constructor cannot find the right core, so we provide
//  A constructur where the user can explicity give the right and
//  left expressions
//
/////////////////////////////////////////////////////////////////////////////

Section::Section(Expression& left_ex_, Expression& right_ex_, int M1, int M2)
  : leftwall_sc(NULL), rightwall_sc(NULL)
{
  global.always_dense = true; // Uniform slabs are not orthogonal here.

  //if ( (M1 < M2) || (M2 < global.N) )
  //  py_print("Warning: M1 > M2 > N not fulfilled.");

  // Add PML and walls to right hand side expression.

  Expression right_ex;
  for (int i=0; i<right_ex_.get_size(); i++)
  {
    Slab* s = dynamic_cast<Slab*>(right_ex_.get_term(i)->get_wg());

    if (s)
    {
      Complex d = right_ex_.get_term(i)->get_d();

      if (i == right_ex_.get_size()-1)
        d += I*global_section.right_PML;

      right_ex += Term((*s)(d));
    }
  }

  if (global_section.rightwall == E_wall)
    rightwall_sc = new E_Wall(*right_ex.get_ext());
  if (global_section.rightwall == H_wall)
    rightwall_sc = new H_Wall(*right_ex.get_ext());
  if (global_section.rightwall != no_wall)
    right_ex.add_term(Term(*rightwall_sc));

  // Add PML and walls to left hand side expression.

  Expression left_ex;
  for (int i=0; i<left_ex_.get_size(); i++)
  {
    Slab* s = dynamic_cast<Slab*>(left_ex_.get_term(i)->get_wg());

    if (s)
    {
      Complex d = left_ex_.get_term(i)->get_d();

      if (i == left_ex_.get_size()-1)
        d += I*global_section.left_PML;

      left_ex += Term((*s)(d));
    }
  }

  if (global_section.leftwall == E_wall)
    leftwall_sc = new E_Wall(*left_ex.get_ext());
  if (global_section.leftwall == H_wall)
    leftwall_sc = new H_Wall(*left_ex.get_ext());
  if (global_section.leftwall != no_wall)
    left_ex.add_term(Term(*leftwall_sc));

  // Create Section.

  s = new Section2D(left_ex, right_ex, M1, M2);
  uniform = s->is_uniform();
  core = s->get_core();

  global.always_dense = false;
}



/////////////////////////////////////////////////////////////////////////////
//
// ModeEstimate::~ModeEstimate
//
/////////////////////////////////////////////////////////////////////////////

ModeEstimate::~ModeEstimate()
{
  delete Ex; Ex = NULL;
  delete Ey; Ey = NULL;
  delete Hx; Hx = NULL;
  delete Hy; Hy = NULL;
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

  sort = highest_index;

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

  Complex z = 0;

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

  Complex z_back = discontinuities.back();
  discontinuities.pop_back();
  discontinuities.push_back(z_back);

  f = new SectionDisp(left,right,global.lambda,M2,symmetric); // TMP 
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
  sort      = section.sort;
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
  sort      = section.sort;
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
// Section2D::material_at
//
/////////////////////////////////////////////////////////////////////////////

Material* Section2D::material_at(const Coord& c) const
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

  return in_left ? left.material_at(c_new) : right.material_at(c_new);
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes()
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

  // If we already calculated modes for a different wavelength/gain
  // combination, use these as an initial estimate, else find them
  // from scratch.

  if (global.sweep_from_previous && (modeset.size() == global.N))
    find_modes_by_sweep();
  else
    find_modes_from_estimates();

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat_n = global.gain_mat->n();
}



/////////////////////////////////////////////////////////////////////////////
//
// Small function objects
//
/////////////////////////////////////////////////////////////////////////////
  
struct index_sorter // Remove?
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
    {
      return ( real(beta_a) > real(beta_b) );
    }
};

struct loss_sorter // Remove?
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
    {
      return ( abs(imag(sqrt(beta_a))) < abs(imag(sqrt(beta_b))) );
    }
};

struct kz2_sorter
{
    bool operator()(const ModeEstimate* a, const ModeEstimate* b)
      {return ( real(a->kz2) > real(b->kz2) );}
};

struct kt_to_neff : ComplexFunction
{
  Complex C;

  kt_to_neff(const Complex& C_) : C(C_) {}

  Complex operator()(const Complex& kt)
    {return sqrt(C - kt*kt)/2./pi*global.lambda;}
};

struct kt2_to_neff : ComplexFunction
{
  Complex C;

  kt2_to_neff(const Complex& C_) : C(C_) {}

  Complex operator()(const Complex& kt2)
    {return sqrt(C - kt2)/2./pi*global.lambda;}
};



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::estimate_kz2_omar_schuenemann
//
/////////////////////////////////////////////////////////////////////////////

vector<ModeEstimate*> Section2D::estimate_kz2_omar_schuenemann()
{
  int n = M1/2;

  // Calc overlap matrices.

  py_print("Creating estimates...");

  Material ref_mat(1.0);
  RefSection ref(ref_mat, get_width(), get_height(), M1);

  cMatrix O_EE(n,n,fortranArray); cMatrix O_MM(n,n,fortranArray);
  cMatrix O_EM(n,n,fortranArray); cMatrix O_zz(n,n,fortranArray);

  ref.calc_overlap_matrices(this, &O_EE, &O_MM, &O_EM, &O_zz);

  // Form auxiliary matrix.

  const Complex omega = 2*pi/global.lambda * c;

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
  if (global.stability == normal)
    kz2 = eigenvalues(E);
  else
    kz2 = eigenvalues_x(E);

  // Return estimates.

  vector<ModeEstimate*> estimates;
  
  for (int i=1; i<=kz2.rows(); i++)
  {
    ModeEstimate* est = new ModeEstimate(kz2(i));
    estimates.push_back(est);
  }

  std::sort(estimates.begin(), estimates.end(), kz2_sorter());

  /*
  for (int i=0; i<estimates.size(); i++)  
    std::cout << i << " " << sqrt(estimates[i]->kz2)/2./pi*global.lambda 
              << std::endl;
  */
  
  return estimates;
}



/////////////////////////////////////////////////////////////////////////////
//
// rotate_slabs()
//
//  Given a sequence of slabs in the x direction, transform it to
//  a sequence of slabs in the y direction.
//
/////////////////////////////////////////////////////////////////////////////

void rotate_slabs(const vector<Slab*>& slabs, const vector<Complex>& disc,
                  vector<Slab*>* slabs_rot,   vector<Complex>* disc_rot)
{
  // Create list of discontinuities in the y-direction.

  disc_rot->clear();
  disc_rot->push_back(0.0);
   for (int i=0; i<slabs.size(); i++)
  {
    vector<Complex> disc_i = slabs[i]->get_discontinuities();
    
    for (int j=0; j<disc_i.size(); j++)
      disc_rot->push_back(disc_i[j]);
  }

  remove_copies(disc_rot, 1e-9);

  std::sort(disc_rot->begin(), disc_rot->end(), RealSorter());

  // Make the rotated slabs.

  for (int i=0; i<disc_rot->size()-1; i++)
  {
    Expression e;
    for (int j=0; j<slabs.size(); j++)
    {
      Material* m = slabs[j]->material_at(Coord((*disc_rot)[i],0.0,0.0,Plus));
      e.add_term((*m)(disc[j+1]-disc[j]));
    }
    
    Slab* slab_i = new Slab(e);
    slabs_rot->push_back(slab_i);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// rebuild_slabs()
//
// Given a stack with a core that hes been split into two pieces, restore
// the stack by removing all artificial discontinuities.
//
// PB: would there be a safer way to achieve this?
//
/////////////////////////////////////////////////////////////////////////////

void rebuild_slabs(const vector<Slab*>& slabs, const vector<Complex>& disc,
                  vector<Slab*>* slabs_rebuilt, vector<Complex>* disc_rebuilt)
{
  vector<Slab*> slabs_rot;
  vector<Complex> disc_x_rot;
  rotate_slabs(slabs, disc, &slabs_rot, &disc_x_rot);
  rotate_slabs(slabs_rot, disc_x_rot, slabs_rebuilt, disc_rebuilt);
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::create_FG_NT
//
//   Noponen/Turunen formulation
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::create_FG_NT(cMatrix* F, cMatrix* G, int M, int N,
                             const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2*pi/global.lambda;
  bool extend = true;

  // Construct data structures for fourier analysis.

  vector<Complex> disc_x(discontinuities);
  disc_x.insert(disc_x.begin(), 0.0);

  vector<vector<Complex> > disc_y, f_eps, f_inv_eps;
  for (int i=0; i<disc_x.size()-1; i++)
  {
    vector<Complex> disc_y_i(slabs[i]->get_discontinuities());
    disc_y_i.insert(disc_y_i.begin(), 0.0);

    vector<Complex> f_eps_i, f_inv_eps_i;
    for (int j=0; j<disc_y_i.size()-1; j++)
    {
      Coord coord(disc_y_i[j],0,0,Plus);
      
      f_eps_i.    push_back(slabs[i]->eps_at(coord)/eps0); 
      f_inv_eps_i.push_back(eps0/slabs[i]->eps_at(coord));
    }

    disc_y.push_back(disc_y_i);
    f_eps.push_back(f_eps_i);
    f_inv_eps.push_back(f_inv_eps_i);
  }

  // Construct fourier matrices.

  cMatrix eps_(4*M+1,4*N+1,fortranArray);
  eps_ = fourier_2D(disc_x, disc_y, f_eps, 2*M, 2*N, extend);

  cMatrix inv_eps_(4*M+1,4*N+1,fortranArray);
  inv_eps_ = fourier_2D(disc_x, disc_y, f_inv_eps, 2*M, 2*N, extend);
  
  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;
  
  cMatrix eps(MN,MN,fortranArray);
  cMatrix inv_eps(MN,MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);
     
      for (int j=-M; j<=M; j++)
        for (int l=-N; l<=N; l++)
        {
          int i2 = (j+M+1) + (l+N)*(2*M+1);

              eps(i1,i2) = eps_(m-j + 2*M+1, n-l + 2*N+1);
          inv_eps(i1,i2) = inv_eps_(m-j + 2*M+1, n-l + 2*N+1);
        }
    }

  // Calculate alpha and beta vector.

  cVector alpha(MN,fortranArray);
  cVector  beta(MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {      
      int i = (m+M+1) + (n+N)*(2*M+1);

      alpha(i) = alpha0 + m*2.*pi/get_width()/2.;
       beta(i) =  beta0 + n*2.*pi/get_height()/2.;
    }

  // Constuct F matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      (*F)(i1,   i2)    =  alpha(i1) * inv_eps(i1,i2) *  beta(i2);
      (*F)(i1,   i2+MN) = -alpha(i1) * inv_eps(i1,i2) * alpha(i2);
      (*F)(i1+MN,i2)    =   beta(i1) * inv_eps(i1,i2) *  beta(i2);
      (*F)(i1+MN,i2+MN) =  -beta(i1) * inv_eps(i1,i2) * alpha(i2);

      if (i1==i2)
      {
        (*F)(i1,   i2+MN) += k0*k0;
        (*F)(i1+MN,i2)    -= k0*k0;
      }
    }

  // Construct G matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {

      (*G)(i1,   i2)    = (i1==i2) ? -alpha(i1)*beta(i2) : 0.0;
      (*G)(i1,   i2+MN) = -k0*k0 * eps(i1,i2);
      (*G)(i1+MN,i2)    =  k0*k0 * eps(i1,i2);
      (*G)(i1+MN,i2+MN) = (i1==i2) ?  alpha(i1)*beta(i2) : 0.0;

      if (i1==i2)
      {
        (*G)(i1,   i2+MN) += alpha(i1)*alpha(i1);
        (*G)(i1+MN,i2)    -=  beta(i2)*beta(i2);
      }
    }

  std::cout << "Created Noponen/Turunen eigenproblem" << std::endl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::create_FG_li
//
//  Li formulation
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::create_FG_li(cMatrix* F, cMatrix* G, int M, int N,
                             const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2*pi/global.lambda;
  bool extend = true;

  // Construct data structures for fourier analysis.

  vector<Complex> disc_x(discontinuities);
  disc_x.insert(disc_x.begin(), 0.0);

  vector<vector<Complex> > disc_y, f_eps;
  for (int i=0; i<disc_x.size()-1; i++)
  {
    vector<Complex> disc_y_i(slabs[i]->get_discontinuities());
    disc_y_i.insert(disc_y_i.begin(), 0.0);

    vector<Complex> f_eps_i, f_inv_eps_i;
    for (int j=0; j<disc_y_i.size()-1; j++)
    {
      Coord coord(disc_y_i[j],0,0,Plus);
      
      f_eps_i.push_back(slabs[i]->eps_at(coord)/eps0);
    }

    disc_y.push_back(disc_y_i);
    f_eps.push_back(f_eps_i);
  }

  vector<Slab*> slabs_rot;
  vector<Complex> disc_x_rot;
  rotate_slabs(slabs, disc_x, &slabs_rot, &disc_x_rot);

  vector<vector<Complex> > disc_y_rot, f_eps_rot, f_inv_eps_rot;
  for (int i=0; i<disc_x_rot.size()-1; i++)
  {
    vector<Complex> disc_y_rot_i(slabs_rot[i]->get_discontinuities());
    disc_y_rot_i.insert(disc_y_rot_i.begin(), 0.0);

    vector<Complex> f_eps_rot_i, f_inv_eps_rot_i;
    for (int j=0; j<disc_y_rot_i.size()-1; j++)
    {
      Coord coord(disc_y_rot_i[j],0,0,Plus);

      f_eps_rot_i.    push_back(slabs_rot[i]->eps_at(coord)/eps0);
      f_inv_eps_rot_i.push_back(eps0/slabs_rot[i]->eps_at(coord));
    }

    disc_y_rot.push_back(disc_y_rot_i);
    f_eps_rot.push_back(f_eps_rot_i);
    f_inv_eps_rot.push_back(f_inv_eps_rot_i);
  }

  // Construct fourier matrices.

  cMatrix eps_(4*M+1,4*N+1,fortranArray);

  eps_ = fourier_2D(disc_x, disc_y, f_eps, 2*M, 2*N, extend);
  
  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;
  
  cMatrix eps(MN,MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);
     
      for (int j=-M; j<=M; j++)
        for (int l=-N; l<=N; l++)
        {
          int i2 = (j+M+1) + (l+N)*(2*M+1);

          eps(i1,i2) = eps_(m-j + 2*M+1, n-l + 2*N+1);
        }
    }

  cMatrix inv_eps_2(MN,MN,fortranArray);
  if (global.stability != SVD)
    inv_eps_2.reference(invert(eps));
  else
    inv_eps_2.reference(invert_svd(eps));

  // Calculate split fourier transforms.

  cMatrix eps_y_x(MN,MN,fortranArray);
  cMatrix eps_x_y(MN,MN,fortranArray);

  cMatrix t(MN,MN,fortranArray);
  t.reference(fourier_2D_split(disc_x_rot,disc_y_rot,f_eps_rot,N,M,extend));
  if (global.stability != SVD)
    eps_y_x.reference(invert(t));
  else
    eps_y_x.reference(invert_svd(t));

  eps_x_y.reference(fourier_2D_split(disc_x_rot,disc_y_rot,f_inv_eps_rot,
                                     N,M,extend));

  // Calculate alpha and beta vector.

  cVector alpha(MN,fortranArray);
  cVector  beta(MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {      
      int i = (m+M+1) + (n+N)*(2*M+1);

      alpha(i) = alpha0 + m*2.*pi/get_width()/2.;
       beta(i) =  beta0 + n*2.*pi/get_height()/2.;
    }

  // Constuct F matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {

      (*F)(i1,   i2)    =  alpha(i1) * inv_eps_2(i1,i2) *  beta(i2);
      (*F)(i1,   i2+MN) = -alpha(i1) * inv_eps_2(i1,i2) * alpha(i2);
      (*F)(i1+MN,i2)    =   beta(i1) * inv_eps_2(i1,i2) *  beta(i2);
      (*F)(i1+MN,i2+MN) =  -beta(i1) * inv_eps_2(i1,i2) * alpha(i2);

      if (i1==i2)
      {
        (*F)(i1,   i2+MN) += k0*k0;
        (*F)(i1+MN,i2)    -= k0*k0;
      }
    }  

  // Construct G matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      (*G)(i1,   i2)    = (i1==i2) ? -alpha(i1)*beta(i2) : 0.0;
      (*G)(i1,   i2+MN) = -k0*k0 * eps_y_x(i1,i2);
      (*G)(i1+MN,i2)    =  k0*k0 * eps_x_y(i1,i2);
      (*G)(i1+MN,i2+MN) = (i1==i2) ?  alpha(i1)*beta(i2) : 0.0;

      if (i1==i2)
      {
        (*G)(i1,   i2+MN) += alpha(i1)*alpha(i1);
        (*G)(i1+MN,i2)    -=  beta(i2)*beta(i2);
      }
    }

  std::cout << "Created Li eigenproblem" << std::endl;

  // Free rotated slabs.

  for (int i=0; i<slabs_rot.size(); i++)
    delete slabs_rot[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::create_FG_li_biaxial
//
//  Li formulation for biaxial media
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::create_FG_li_biaxial(cMatrix* F, cMatrix* G, int M, int N,
                             const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2*pi/global.lambda;
  bool extend = true;

  const Real p = global_section.PML_fraction;
  Complex p_d = 0.3; // TMP

  // Construct data structures for fourier analysis.

  vector<Complex> disc_x(discontinuities);
  disc_x.insert(disc_x.begin(), 0.0);

  for (int i=0; i<disc_x.size(); i++)
    disc_x[i] = real(disc_x[i]);

  vector<vector<Complex> > disc_y;
  vector<vector<Complex> > f_eps_1, f_eps_2, f_eps_3;
  vector<vector<Complex> > f_mu_1,  f_mu_2,  f_mu_3;

  for (int i=0; i<disc_x.size()-1; i++)
  {
    vector<Complex> disc_y_i(slabs[i]->get_discontinuities());
    disc_y_i.insert(disc_y_i.begin(), 0.0);

    for (int j=0; j<disc_y_i.size(); j++)
      disc_y_i[j] = real(disc_y_i[j]);
    
    vector<Complex> f_eps_1_i, f_eps_2_i, f_eps_3_i;
    vector<Complex>  f_mu_1_i,  f_mu_2_i,  f_mu_3_i;
   
    for (int j=0; j<disc_y_i.size()-1; j++)
    {
      Material* m = slabs[i]->material_at(Coord(disc_y_i[j],0,0,Plus));
      
      f_eps_1_i.push_back(1./m->epsr(1));      
      f_eps_2_i.push_back(   m->epsr(2));
      f_eps_3_i.push_back(   m->epsr(3));      

      f_mu_1_i. push_back(1./m->mur(1));      
      f_mu_2_i. push_back(   m->mur(2));
      f_mu_3_i. push_back(   m->mur(3));
    }

    // Add lower PML

    if (abs(global_slab.lower_PML) > 1e-12)
    {
      Complex d = disc_y_i[1];
      disc_y_i.insert(disc_y_i.begin()+1, p_d);

      Complex a = 1. + global_slab.lower_PML / p_d * I;
      
      f_eps_1_i.insert(f_eps_1_i.begin(), f_eps_1_i[0]/a);
      f_eps_2_i.insert(f_eps_2_i.begin(), f_eps_2_i[0]/a);   
      f_eps_3_i.insert(f_eps_3_i.begin(), f_eps_3_i[0]*a);

      f_mu_1_i .insert(f_mu_1_i. begin(), f_mu_1_i[0]/a);    
      f_mu_2_i .insert(f_mu_2_i. begin(), f_mu_2_i[0]/a);   
      f_mu_3_i .insert(f_mu_3_i. begin(), f_mu_3_i[0]*a);
    }
    
    // Add upper PML.

    if (abs(global_slab.upper_PML) > 1e-12)
    {      
      Complex d = disc_y_i[disc_y_i.size()-1] - disc_y_i[disc_y_i.size()-2];
      disc_y_i.insert(disc_y_i.end()-1, disc_y_i.back() - p_d);

      Complex a = 1. + global_slab.upper_PML / (p_d) * I;

      f_eps_1_i.push_back(f_eps_1_i.back()/a);    
      f_eps_2_i.push_back(f_eps_2_i.back()/a);   
      f_eps_3_i.push_back(f_eps_3_i.back()*a);

      f_mu_1_i .push_back(f_mu_1_i.back()/a);    
      f_mu_2_i .push_back(f_mu_2_i.back()/a);
      f_mu_3_i .push_back(f_mu_3_i.back()*a);
    }

    // Finalise.

    disc_y.push_back(disc_y_i);

    f_eps_1.push_back(f_eps_1_i);    
    f_eps_2.push_back(f_eps_2_i);   
    f_eps_3.push_back(f_eps_3_i);

    f_mu_1 .push_back(f_mu_1_i);    
    f_mu_2 .push_back(f_mu_2_i);   
    f_mu_3 .push_back(f_mu_3_i);
  }

  // Add left PMl.

  if (abs(global_section.left_PML) > 1e-12)
  {
    Complex d = disc_x[1];
    disc_x.insert(disc_x.begin()+1, p_d);

    vector<Complex> f_eps_1_i, f_eps_2_i, f_eps_3_i;
    vector<Complex>  f_mu_1_i,  f_mu_2_i,  f_mu_3_i;

    Complex a = 1. + global_section.left_PML / p_d * I;

    for (int j=0; j<disc_y[0].size()-1; j++)
    {
      f_eps_1_i.push_back(f_eps_1[0][j]*a);
      f_eps_2_i.push_back(f_eps_2[0][j]*a);
      f_eps_3_i.push_back(f_eps_3[0][j]*a);
      
      f_mu_1_i .push_back( f_mu_1[0][j]*a);
      f_mu_2_i .push_back( f_mu_2[0][j]*a);
      f_mu_3_i .push_back( f_mu_3[0][j]*a);
    }
    
    disc_y.insert(disc_y.begin(), disc_y[0]);

    f_eps_1.insert(f_eps_1.begin(), f_eps_1_i);
    f_eps_2.insert(f_eps_2.begin(), f_eps_2_i);
    f_eps_3.insert(f_eps_3.begin(), f_eps_3_i);
  
    f_mu_1 .insert(f_mu_1 .begin(), f_mu_1_i);    
    f_mu_2 .insert(f_mu_2 .begin(), f_mu_2_i);   
    f_mu_3 .insert(f_mu_3 .begin(), f_mu_3_i);
  }

  // Add right PMl.

  if (abs(global_section.right_PML) > 1e-12)
  {
    Complex d = disc_x[disc_x.size()-1] - disc_x[disc_x.size()-2];
    disc_x.insert(disc_x.end()-1, disc_x.back() - p_d);

    vector<Complex> f_eps_1_i, f_eps_2_i, f_eps_3_i;
    vector<Complex>  f_mu_1_i,  f_mu_2_i,  f_mu_3_i;

    Complex a = 1. + global_section.right_PML / p_d * I;

    for (int j=0; j<disc_y.back().size()-1; j++)
    {
      f_eps_1_i.push_back(f_eps_1.back()[j]*a);
      f_eps_2_i.push_back(f_eps_2.back()[j]*a);
      f_eps_3_i.push_back(f_eps_3.back()[j]*a);
      
      f_mu_1_i .push_back( f_mu_1.back()[j]*a);
      f_mu_2_i .push_back( f_mu_2.back()[j]*a);
      f_mu_3_i .push_back( f_mu_3.back()[j]*a);
    }
    
    disc_y.push_back(disc_y.back());

    f_eps_1.push_back(f_eps_1_i);
    f_eps_2.push_back(f_eps_2_i);
    f_eps_3.push_back(f_eps_3_i);
  
    f_mu_1 .push_back( f_mu_1_i);    
    f_mu_2 .push_back( f_mu_2_i);   
    f_mu_3 .push_back( f_mu_3_i);
  }

  // Construct fourier matrices.

  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;

  cMatrix eps_1(MN,MN,fortranArray); cMatrix eps_2_(MN,MN,fortranArray);
  cMatrix  mu_1(MN,MN,fortranArray); cMatrix  mu_2_(MN,MN,fortranArray);

  eps_1. reference(fourier_2D_split(disc_x, disc_y, f_eps_1, M, N, extend));
   mu_1. reference(fourier_2D_split(disc_x, disc_y,  f_mu_1, M, N, extend)); 
 
  eps_2_.reference(fourier_2D_split(disc_x, disc_y, f_eps_2, M, N, extend));
   mu_2_.reference(fourier_2D_split(disc_x, disc_y,  f_mu_2, M, N, extend));  

  cMatrix eps_3_f(4*M+1,4*N+1,fortranArray); 
  cMatrix  mu_3_f(4*M+1,4*N+1,fortranArray);
  
  eps_3_f = fourier_2D(disc_x, disc_y, f_eps_3, 2*M, 2*N, extend);
   mu_3_f = fourier_2D(disc_x, disc_y,  f_mu_3, 2*M, 2*N, extend);

  cMatrix eps_3_(MN,MN,fortranArray), mu_3_(MN,MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (n+N+1) + (m+M)*(2*N+1);
     
      for (int j=-M; j<=M; j++)
        for (int l=-N; l<=N; l++)
        {
          int i2 = (l+N+1) + (j+M)*(2*N+1);

          eps_3_(i1,i2) = eps_3_f(m-j + 2*M+1, n-l + 2*N+1);          
           mu_3_(i1,i2) =  mu_3_f(m-j + 2*M+1, n-l + 2*N+1);
        }
    }

  cMatrix eps_2(MN,MN,fortranArray), mu_2(MN,MN,fortranArray); 
  cMatrix eps_3(MN,MN,fortranArray), mu_3(MN,MN,fortranArray);

  if (global.stability != SVD)
  {
    eps_2.reference(invert(eps_2_));    
    eps_3.reference(invert(eps_3_));

    mu_2. reference(invert(mu_2_));
    mu_3. reference(invert(mu_3_));
  }
  else
  {
    eps_2.reference(invert_svd(eps_2_));    
    eps_3.reference(invert_svd(eps_3_));

    mu_2. reference(invert_svd(mu_2_));
    mu_3. reference(invert_svd(mu_3_));
  }

  // Calculate alpha and beta vector.

  cVector alpha(MN,fortranArray);
  cVector  beta(MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {      
      int i = (m+M+1) + (n+N)*(2*M+1);

      alpha(i) = alpha0 + m*2.*pi/real(get_width())/2.;
       beta(i) =  beta0 + n*2.*pi/real(get_height())/2.;
    }

  // Constuct F matrix.  

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      (*F)(i1,   i2)    =  alpha(i1) * eps_3(i1,i2) *  beta(i2);
      (*F)(i1,   i2+MN) = -alpha(i1) * eps_3(i1,i2) * alpha(i2);
      (*F)(i1+MN,i2)    =   beta(i1) * eps_3(i1,i2) *  beta(i2);
      (*F)(i1+MN,i2+MN) =  -beta(i1) * eps_3(i1,i2) * alpha(i2);

      (*F)(i1,   i2+MN) += k0*k0 * mu_2(i1,i2);
      (*F)(i1+MN,i2)    -= k0*k0 * mu_1(i1,i2);
    }

  // Construct G matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      (*G)(i1,   i2)    = -alpha(i1) * mu_3(i1,i2) *  beta(i2);
      (*G)(i1,   i2+MN) =  alpha(i1) * mu_3(i1,i2) * alpha(i2);
      (*G)(i1+MN,i2)    =  -beta(i1) * mu_3(i1,i2) *  beta(i2);
      (*G)(i1+MN,i2+MN) =   beta(i1) * mu_3(i1,i2) * alpha(i2);

      (*G)(i1,   i2+MN) -= k0*k0 * eps_2(i1,i2);
      (*G)(i1+MN,i2)    += k0*k0 * eps_1(i1,i2);
    }

  std::cout << "Created Li biaxial eigenproblem" << std::endl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::create_FG_ASR
//
// Adaptive Spatial Resolution formulation of the Li algorithm.
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::create_FG_ASR(cMatrix* F, cMatrix* G, int M, int N,
                              const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2*pi/global.lambda;
  bool  extend     = true;
  const Real eta   = global_section.eta_ASR;
  bool  stretching = (global_section.section_solver == ASR_2D_stretched);
  
  // These booleans change the integration sequence.

  bool A_switch = global_section.A_switch;
  bool B_switch = global_section.B_switch;
  bool C_switch = global_section.C_switch;
  bool D_switch = global_section.D_switch;
  
  // Varying these parameters has shown that only the combination 
  // false, true, false, true works.

  // [PDB] this is kind of arbitrary here, I just wanted to see the effect 
  // of changing the stretching parameters.

  const Real u_step_given = global_section.u_step_given;
  const Real v_step_given = global_section.v_step_given;

  /*
  [PDB]
  Again, this is a desperate measure in order to find the error in 
  the algorithm, this value is a value between 1 and 0
  ---> 0 no stretching, so only ASR
  ---> 1 full stretching, ASR_stretched
  --->everything in between is, well, in between...
  This method is only to be used in combination with asr_2D_stretched
  */
  
  const Real percentage_stretched = global_section.percentage_stretched;
  bool extended_output = global_section.extended_output;

  // Construct data structures for fourier analysis.

  // [PDB] Because there is a problem with the artificial discontinuity
  // in the core, we remove this discontinuity with the function rebuild_slabs.
  
  // Preliminary declaration of disc_x (with artificial core discontinuity).
  
  vector<Complex> disc_x_temp(discontinuities);
  disc_x_temp.insert(disc_x_temp.begin(), 0.0);
  
  // Construct slab structure without an artificial discontinuity in the core.

  vector<Slab*> slabs_new;
  vector<Complex> disc_x;

  rebuild_slabs(slabs, disc_x_temp, &slabs_new, &disc_x);

  // Create list of discontinuities in the y direction.

  // Search for largest possible section in the y direction (slab-thickness)
  // PB: I don't understand this comment. disc_y_max --> disc_y?
  
  vector<Complex> disc_y_max;
  disc_y_max.push_back(0.0);
  
  for (int i=0; i<slabs_new.size(); i++)
  {
    vector<Complex> disc_y_i = slabs_new[i]->get_discontinuities();
    
    for (int j=0; j<disc_y_i.size(); j++)
      disc_y_max.push_back(disc_y_i[j]);
  }

  remove_copies(&disc_y_max, 1e-9);
  std::sort(disc_y_max.begin(), disc_y_max.end(), RealSorter());
 
  // Create disc_x, disc_y, f_eps, f_inv_eps.	.
  // While doing this, look for step_u and step_v in order to save some time.

  Real step_u = 0.0;
  Real step_v = 0.0;

  if (extended_output)
    for (int i=0; i<disc_x.size(); i++)
      std::cout << "disc_x(" << i << ") = " << disc_x[i] << std::endl;

  // Create extra matrix h_eps, which has the same dimension as f_eps but 
  // which is filled with the value 1.0.

  vector<vector<Complex> > disc_y, f_eps, f_inv_eps, h_eps;
  for (int i=0; i<disc_x.size()-1; i++)
  {
     // Determine biggest step in x direction for ASR stretching algorithm.

    if (real(disc_x[i+1]-disc_x[i]) > step_u)
      step_u = real(disc_x[i+1] - disc_x[i]);

    //vector<Complex> disc_y_i(slabs[i]->get_discontinuities());
    //disc_y_i.insert(disc_y_i.begin(), 0.0);
    // PB: Why can't you use this?
   
    vector<Complex> f_eps_i,f_inv_eps_i,h_eps_i;
    for (int j=0; j<disc_y_max.size()-1; j++)
    {
       // Determine biggest step in y direction for ASR stretching algorithm.

      if (real(disc_y_max[j+1] - disc_y_max[j])> step_v)
        step_v = real(disc_y_max[j+1] - disc_y_max[j]);

      Coord coord(disc_y_max[j],0,0,Plus);
      
      f_eps_i.push_back(slabs_new[i]->eps_at(coord)/eps0);
      f_inv_eps_i.push_back(eps0/slabs_new[i]->eps_at(coord));
      h_eps_i.push_back(1.0); 	
    }

    disc_y.push_back(disc_y_max); // PB: all columns same?
    f_eps.push_back(f_eps_i);
    f_inv_eps.push_back(f_inv_eps_i);
    h_eps.push_back(h_eps_i);
  }

  if (extended_output)
  {
    for (int i=0; i<disc_x.size()-1; i++)
      for (int j=0; j<disc_y[i].size(); j++)
        std::cout << "disc_y(" << i << "," << j << ") = " 
                  << disc_y[i][j] <<std::endl;

    for (int i=0; i<disc_x.size()-1; i++)
      for (int j=0; j<f_eps[i].size(); j++)
        std::cout << "f_eps(" << i << "," << j << ") = " 
                  << f_eps[i][j] << std::endl;

    for (int i=0; i<disc_x.size()-1; i++)
      for (int j=0; j<f_inv_eps[i].size(); j++)
        std::cout << "f_inv_eps(" << i << "," << j <<" ) = " 
                  << f_inv_eps[i][j] <<std::endl;

    for (int i=0; i<disc_x.size()-1; i++)
      for (int j=0; j<h_eps[i].size(); j++)
        std::cout << "h_eps(" << i << "," << j << ") = "
                  << h_eps[i][j] << std::endl;
  }

  // Create rotated slabs.

  vector<Slab*> slabs_rot;
  vector<Complex> disc_x_rot;
  rotate_slabs(slabs_new, disc_x, &slabs_rot, &disc_x_rot);
  
  // Create list of discontinuities in the y-direction.

  vector<Complex> disc_y_rot_max;
  disc_y_rot_max.clear();
  disc_y_rot_max.push_back(0.0);
  
  for (int i=0; i<slabs_rot.size(); i++)
  {
    vector<Complex> disc_y_i = slabs_rot[i]->get_discontinuities();
    
    for (int j=0; j<disc_y_i.size(); j++)
      disc_y_rot_max.push_back(disc_y_i[j]);
  }

  remove_copies(&disc_y_rot_max, 1e-9);
  std::sort(disc_y_rot_max.begin(), disc_y_rot_max.end(), RealSorter());
 
  // Again we have to create a matrix with the same dimensions as f_eps_rot 
  // but filled with the value 1.0.

  vector<vector<Complex> > disc_y_rot, f_eps_rot, f_inv_eps_rot,h_eps_rot;
  for (int i=0; i<disc_x_rot.size()-1; i++)
  { 
    //vector<Complex> disc_y_rot_i(slabs_rot[i]->get_discontinuities());
    //disc_y_rot_i.insert(disc_y_rot_i.begin(), 0.0);

    vector<Complex> f_eps_rot_i,f_inv_eps_rot_i,h_eps_rot_i;
    for (int j=0; j<disc_y_rot_max.size()-1; j++)
    {
      Coord coord(disc_y_rot_max[j],0,0,Plus);

      // PDB: Changed this statement.

      f_eps_rot_i.push_back(slabs_rot[i]->eps_at(coord)/eps0);
      f_inv_eps_rot_i.push_back(eps0/slabs_rot[i]->eps_at(coord));
      h_eps_rot_i.push_back(1.0);
    }

    disc_y_rot.push_back(disc_y_rot_max);
    f_eps_rot.push_back(f_eps_rot_i);
    f_inv_eps_rot.push_back(f_inv_eps_rot_i);
    h_eps_rot.push_back(h_eps_rot_i);
  }

  if (extended_output)
  {
    for(int i=0; i<disc_x_rot.size()-1; i++)
      for(int j=0; j<disc_y_rot[i].size(); j++)
        std::cout << "disc_y_rot(" << i << "," << j << ") = "
                  << disc_y_rot[i][j] << std::endl;

    for(int i=0; i<disc_x_rot.size()-1; i++)
      for(int j=0; j<f_inv_eps_rot[i].size(); j++)
        std::cout << "f_inv_eps_rot(" << i << "," << j << ") = "
                  << f_inv_eps_rot[i][j] << std::endl;

    for(int i=0; i<disc_x_rot.size()-1; i++)
      for(int j=0; j<f_eps_rot[i].size(); j++)
        std::cout << "f_eps_rot(" << i << "," << j << ") = "
                  << f_eps_rot[i][j] << std::endl;

    for(int i=0; i<disc_x_rot.size()-1; i++)
      for(int j=0; j<h_eps_rot[i].size(); j++)
        std::cout << "h_eps_rot(" << i << "," << j << ") = "
                  << h_eps_rot[i][j] << std::endl;	
  }

  /*

  CODE EXPLANATION

  Construct the ASR versions of the matrices 
  x_disc, y_disc, x_disc_rot and y_disc_rot.

  Method:	 

  1. Find the biggest thickness of the structure along x
  2. Create a number of equally thick sections along x, with the thickness 
     of the thickest structure
  3. Find the thickest structure along the y direction
  4. Create a number of equally thick sections along y, with the thickness 
     of the thickest structure
  
  Step 1 and 3 have already been performed earlier in order to save time.

  The following code is actually just an experiment
  How about taking step_u = step_v = the average of both
  Maybe this will make the ASR algorithm more consistent

  Real step = (step_u + step_v)/2.;
  std::cout << "Stepsizes " << " u " << step_u << " v " << step_v << " avg "
            << step << std::endl;
  step_u = step;
  step_v = step;

  I've tried this but there is no improvement compared to the old formulation

  */

  // [PDB] The following code lets you chose the amount of stretching by
  // chosing step_u and step_v yourself.

  if (abs(u_step_given) > 1e-16)
  {
   step_u = u_step_given;
   std::cout << "Using the given step for u" <<std::endl;
  }
 
  if (abs(v_step_given) > 1e-16)
  {
   step_v = v_step_given;
   std::cout << "Using the given step for v" << std::endl;
  }
 
  // Creation of the new disc_u and disc_v vector.
  // I also took care of the PML conditions (I think).

  vector<Complex> disc_u;
  disc_u.insert(disc_u.begin(), 0.0);
 
  // Question: why doesn't the first element have a PML contribution.

  vector<vector<Complex> > disc_v;

  // Checking the ASR_stretched algorithm in general.

  if (stretching)
  {
    for (int i=1; i<disc_x.size(); i++)
    {
      disc_u.push_back(percentage_stretched*Real(i)*step_u 
                       +(1.-percentage_stretched)*
                       real(disc_x[i])+ imag(disc_x[i])*I);

      // [PDB] I have messed with this code in order to be able to change 
      // the stretching percentage.

      vector<Complex> disc_v_i;
      disc_v_i.insert(disc_v_i.begin(), 0.0);
      for (int j=1; j<disc_y[i-1].size(); j++)
        disc_v_i.push_back(percentage_stretched*Real(j)*step_v
                           +(1.-percentage_stretched)*real(disc_y[i-1][j]) 
                           + imag(disc_y[i-1][j])*I);
      
      disc_v.push_back(disc_v_i);
    }
  }
  else
  {
    disc_u = disc_x;
    disc_v = disc_y;
  }

  // Question: Why doesn't the first element have a PML contribution?

  if (extended_output)
  {
    for (int i=0; i<disc_u.size(); i++)
      std::cout << "disc_u(" << i << ") = " << disc_u[i] << std::endl;

    for (int i=0; i<disc_u.size()-1; i++)
      for (int j=0; j<disc_v[i].size(); j++)
        std::cout << "disc_v(" << i << "," << j << ") = "
                  << disc_v[i][j] << std::endl;
  }

  // Create rotated matrices disc_u_rot, disc_v_rot.

  vector<Complex> disc_u_rot;
  disc_u_rot.insert(disc_u_rot.begin(), 0.0);

  // Question: why doesn't the first element have a PML contribution?

  vector<vector<Complex> > disc_v_rot;
  if (stretching)
  {
    for (int i=1; i<(disc_v[0]).size(); i++)
    {
      disc_u_rot.push_back(percentage_stretched*Real(i)*step_v 
                           + (1.-percentage_stretched)*real(disc_y[0][i]) 
                           + imag(disc_v[0][i])*I);

      vector<Complex> disc_v_rot_i;
      disc_v_rot_i.insert(disc_v_rot_i.begin(), 0.0);
      for (int j=1; j<disc_u.size(); j++)
        disc_v_rot_i.push_back(percentage_stretched*Real(j)*step_u 
                               + (1.-percentage_stretched)*real(disc_x[j])
                               + imag(disc_x[j])*I);

      disc_v_rot.push_back(disc_v_rot_i);
    }
  }
  else
  {
    disc_u_rot = disc_x_rot;
    disc_v_rot = disc_y_rot;
  }

  if (extended_output)
  {
    for (int i=0; i<disc_u_rot.size(); i++)
      std::cout << "disc_u_rot("<<i<<") = " << disc_u_rot[i] << std::endl;

    for (int i=0; i<disc_u_rot.size()-1; i++)
      for (int j=0; j<disc_v_rot[i].size(); j++)
        std::cout << "disc_v_rot(" << i << "," << j << ") = " 
                  << disc_v_rot[i][j] << std::endl;
  }

  /*
  We have now constructed the following data structures
  TODO_PDB exploit the fact that some matrices only have to be calculated 
  for Li_1 or for Li_2
  
  disc_x        Discontinuities along the x-axis
  disc_y        Discontinuities along the y-axis 
  f_eps	        Dielectric constant
  f_inv_eps     Inverse of the dielectric constant, ONLY FOR LI_1
  h_eps	        Matrix with same dimension as f_eps, filed with 1.
  
  disc_x_rot    Discontinuities along the x-axis, rotated structure
  disc_y_rot    Discontinuities along the y-axis, rotated sytructure
  f_eps_rot     Dielectric constant, rotated structure, ONLY FOR LI_2
  f_inv_eps_rot	Inverse of the dielectric constant, rotated structure
  h_eps_rot     Matrix with same dimension as f_inv_eps_rot, filled with 1

  disc_u	Discontinuities along the u-axis
  disc_v	Discontinuities along the v-axis
  
  disc_u_rot	Discontinuities along the u-axis, rotated structure
  disc_v_rot    Discontinuities along the v-axis, rotated structure
  */

  // Construct fourier matrices
	
  // [PDB] First we will calculate the Fourier matrices for which the Li 
  // formulation does not apply.

  // Normal matrices.
  
  cMatrix eps_f_g_(4*M+1,4*N+1,fortranArray);
  cMatrix     f_g_(4*M+1,4*N+1,fortranArray);

  eps_f_g_ = fourier_2D_ASR(disc_x,disc_y,disc_u,disc_v,f_eps,
                            2*M,2*N,extend,eta);
  f_g_	   = fourier_2D_ASR(disc_x,disc_y,disc_u,disc_v,h_eps,
                            2*M,2*N,extend,eta);
  
  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;
  
  // (Double) Toeplitz matrices.

  cMatrix eps_f_g(MN,MN,fortranArray);
  cMatrix f_g(MN,MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);
     
      for (int j=-M; j<=M; j++)
        for (int l=-N; l<=N; l++)
        {
          int i2 = (j+M+1) + (l+N)*(2*M+1);

          eps_f_g(i1,i2) = eps_f_g_(m-j + 2*M+1, n-l + 2*N+1);
          f_g(i1,i2) = f_g_(m-j + 2*M+1, n-l + 2*N+1);
        }
    }

  // Inverse (double) Toeplitz matrices.

  cMatrix inv_eps_f_g(MN,MN,fortranArray);
  cMatrix inv_f_g    (MN,MN,fortranArray);
  
  if (global.stability != SVD)
  {
    inv_eps_f_g.reference(invert(eps_f_g));
    inv_f_g.    reference(invert(f_g));
  }   
  else
  {
    inv_eps_f_g.reference(invert_svd(eps_f_g));
    inv_f_g.    reference(invert_svd(f_g));
  }
  
  if (extended_output)
  {
    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "inv_eps_f_g(" << i << "," << j << ") = "
                  << inv_eps_f_g(i,j) <<std::endl;

    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "inv_f_g(" << i << "," << j << ") = "
                  << inv_f_g(i,j) <<std::endl;
  }

  /*
  Now we will calculate the Fourier matrices for which the Li formulation 
  does apply. We wil use fourier_2D_split_ASR instead of fourier_2D_ASR
  
  Note:
  Li uses a somewhat different form for the Fourier matrices in the paper
  "Fourier modal methods for crossed anisotropic gratings with arbitrary 
  permittivity and permeability tensors", Josa A, 5(2002), 345-355	
  However, the small difference in convergence is a purely numerical effect, 
  the faster convergence of this method for the ASR configuration remains 
  unclear up till now. That is why I have included all possible formulations 
  in this solver. It is possible to switch  between formulations by changing 
  the booleans A_switch, B_switch, C_switch and D_switch.
 
  These booleans change the integration sequence of the matrices
  g_f_v_u                               
  inv_eps_g_f_v_u
  f_g_u_v
  inv_eps_f_g_u_v

  */

  cMatrix         g_f_v_u(MN,MN,fortranArray);
  cMatrix inv_eps_g_f_v_u(MN,MN,fortranArray);
  cMatrix         f_g_u_v(MN,MN,fortranArray);
  cMatrix inv_eps_f_g_u_v(MN,MN,fortranArray);
   
  // TODO_PDB Simplify fourier calculations of gfvu and fguv, by not using 
  // pseudo fourier transforms but full Fourier transforms in both directions.
  
  // Calculation of the matrix F_G_U_V, so what I always refer to as A

  if (!A_switch)
  {
    // Calculation of the fourier components of f according to u, and of g
    // according to v, this routine uses the rotated structure.

    f_g_u_v.reference(fourier_2D_split_ASR(disc_x_rot,disc_y_rot,disc_u_rot,
                      disc_v_rot,h_eps_rot,N,M,extend,eta));
  }
  else
  {
    // For simplicity reasons we will keep the names that were assigned to 
    // these matrice for the original Li formulation, although it is quite 
    // obvious now that these names do not reflect the correct mathematical 
    // approach for the modified method.

    // Calculation of the fourier components of g according to v, and of f
    // according to u, this routine uses the regular structure.

    cMatrix temp(MN,MN,fortranArray);

    temp.reference(fourier_2D_split_ASR(disc_x,disc_y,disc_u,disc_v,
                                        h_eps,N,M,extend,eta));

    if (global.stability != SVD)
      f_g_u_v.reference(invert(temp));
    else
      f_g_u_v.reference(invert_svd(temp));
  }

  // Calculation of the matrix G_F_V_U, so what I always refer to as B.

  if (!B_switch) // Original.
  {
    // Calculation of the fourier components of g according to v, 
    // and of f according to u.

    g_f_v_u.reference(fourier_2D_split_ASR(disc_x,disc_y,disc_u,disc_v,
                                           h_eps,M,N,extend,eta));
  }
  else
  {
    // Calculation of the fourier components of f according to u and of g 
    // according to v.

    cMatrix temp(MN,MN,fortranArray);

    temp.reference(fourier_2D_split_ASR(disc_x_rot,disc_y_rot,disc_u_rot, 
                                        disc_v_rot,h_eps_rot,N,M,extend,eta));

    if (global.stability != SVD)
      g_f_v_u.reference(invert(temp));
    else
      inv_eps_g_f_v_u.reference(invert_svd(temp));
  }

  // Calculation of the matrix inv_eps_f_g_u_v, so what I always refer to as C

  if (!C_switch)
  {
    inv_eps_f_g_u_v.reference(fourier_2D_split_ASR(disc_x_rot,disc_y_rot,
                                                   disc_u_rot,disc_v_rot,
                                                   f_inv_eps_rot,N,M,extend, 
                                                   eta));
  }
  else
  {
    // Calculation of the fourier components of g/eps according to v, and of f
    // according to u, this routine uses the regular structure.

    cMatrix temp(MN,MN,fortranArray);

    temp.reference(fourier_2D_split_ASR(disc_x,disc_y,disc_u,disc_v,
                                        f_inv_eps,N,M,extend,eta));

    if (global.stability != SVD)
      inv_eps_f_g_u_v.reference(invert(temp));
    else
      inv_eps_f_g_u_v.reference(invert_svd(temp));
   }

  // Calculation of the matrix INV_EPS_G_F_V_U, so what I always refer to as D
  
  if (!D_switch) 
  {
    // Calculation of the fourier components of g/eps according to v, and of 
    // f according to u.

    inv_eps_g_f_v_u.reference(fourier_2D_split_ASR(disc_x,disc_y,disc_u,
                                                   disc_v,f_inv_eps,
                                                   M,N,extend,eta));
  }
  else
  {
    // Calculation of the fourier components of eps*f according to u, and of g
    // according to v.

    cMatrix temp(MN,MN,fortranArray);

    temp.reference(fourier_2D_split_ASR(disc_x_rot,disc_y_rot,
                                        disc_u_rot,disc_v_rot,
                                        f_eps_rot,N,M,extend,eta));

    if (global.stability != SVD)
      inv_eps_g_f_v_u.reference(invert(temp));
    else
      inv_eps_g_f_v_u.reference(invert_svd(temp));
  }

  if (extended_output)
  {
    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "g_f_v_u(" << i << "," << j << ") = "
                  << g_f_v_u(i,j) << std::endl;

    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "inv_eps_g_f_v_u(" << i << "," << j << ") ="
                  << inv_eps_g_f_v_u(i,j) << std::endl;

    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "f_g_u_v(" << i << "," << j << ") =" 
                  << f_g_u_v(i,j) << std::endl;

    for (int i=1; i<=MN; i++)
      for (int j=1; j<=MN; j++)
        std::cout << "inv_eps_f_g_u_v(" << i << "," << j << ")  ="
                  << inv_eps_f_g_u_v(i,j) << std::endl;
  }

  // Calculate alpha and beta vector.

  cVector alpha(MN,fortranArray);
  cVector  beta(MN,fortranArray);

  /*
  Since the structure will be stretched(sometimes), the functions get_width()
  and get_height() cannot be trusted anymore. Therefore we calculate here the
  exact (stretched) height and width
  */

  const Complex Lu = disc_u.back()     - disc_u.front();
  const Complex Lv = disc_u_rot.back() - disc_u_rot.front();

  if (extended_output)
    std::cout << "Lu " << Lu << ", Lv " << Lv << std::endl;

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {      
      int i = (m+M+1) + (n+N)*(2*M+1);

      alpha(i) = alpha0 + m*2.*pi/Lu/2.; // This used to be get_width()
       beta(i) =  beta0 + n*2.*pi/Lv/2.; // This used to be get_height()
    }
  
  if (extended_output)
  {
    for (int i=1; i<=MN; i++)
      std::cout << "alpha(" << i << ") = " << alpha(i) << std::endl;

    for (int i=1; i<=MN; i++)
      std::cout << "beta(" << i << ") = " << beta(i) << std::endl;
  }

  // Construction of the F and G matrix of the ASR formalism
  // Construct F matrix.

  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {

      (*F)(i1,   i2)    =  alpha(i1) * inv_eps_f_g(i1,i2) *  beta(i2);
      (*F)(i1,   i2+MN) =  k0 * k0 * g_f_v_u(i1,i2) 
                            - alpha(i1) * inv_eps_f_g(i1,i2) * alpha(i2);
      (*F)(i1+MN,i2)    = -k0 * k0 * f_g_u_v(i1,i2) 
                            + beta(i1) *inv_eps_f_g(i1,i2) * beta(i2);
      (*F)(i1+MN,i2+MN) = -beta(i1) * inv_eps_f_g(i1,i2) * alpha(i2);

    }

  // Construct G matrix.
	
  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      //if (eta == 0.0)
      //  (*G)(i1,   i2)    = (i1==i2) ? -alpha(i1)*beta(i2) : 0.0;
      //else

	  (*G)(i1,   i2)    = -alpha(i1) * inv_f_g(i1,i2) * beta(i2);
	  (*G)(i1,   i2+MN) = -k0 * k0 * inv_eps_g_f_v_u(i1,i2) 
                                + alpha(i1) * inv_f_g(i1,i2) * alpha(i2);
          (*G)(i1+MN,i2)    =  k0 * k0 * inv_eps_f_g_u_v(i1,i2) 
                                - beta(i1) * inv_f_g(i1,i2) * beta(i2);

      //if (eta == 0.0)
      //  (*G)(i1+MN,i2+MN) = (i1==i2) ?  alpha(i1)*beta(i2) : 0.0;
      //else
          (*G)(i1+MN,i2+MN) =  alpha(i1) * inv_f_g(i1,i2) * beta(i2) ;  
     
    }

    /*

    In case of eta=0 the ASR formalism should be reduced to the Li formalism, 
    that is why we check here what happend to the matrices inv_f_g
    
    for (int i=1;i<MN;i++)
      for (int j=1;j<MN;j++)
        std::cout << "inv_f_g(" << i << "," << j << ") =" 
                  << real(inv_f_g(i,j)) << " " << imag(inv_f_g(i,j))
                  << std::endl;
   */


  if (extended_output)
  {
    for (int i=1; i<=2*MN; i++)
      for (int j=1; j<=2*MN; j++)
        std::cout << "F(" << i << "," << j <<") ="
                  << real((*F)(i,j)) << " " << imag((*F)(i,j)) <<std::endl;

    for (int i=1; i<=2*MN; i++)
      for (int j=1; j<=2*MN; j++)
        std::cout << "G(" << i << "," << j <<") ="
                  << real((*G)(i,j)) << " "<< imag((*G)(i,j)) << std::endl;
  }

  std::cout << "Created ASR eigenproblem" << std::endl;

  // Free rotated slabs.
  
  for (int i=0; i<slabs_rot.size(); i++)
    delete slabs_rot[i];

  // Free rebuilt slabs.
  
  for (int i=0; i<slabs_new.size(); i++)
    delete slabs_new[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// IndexConvertor
//
/////////////////////////////////////////////////////////////////////////////

struct IndexConvertor
{
    IndexConvertor(int M_, int N_): M(M_), N(N_) {}
    
    int operator() (int m, int n) const {return (m+M+1) + (n+N)*(2*M+1);}

    int M,N;
};



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::estimate_kz2_fourier
//
/////////////////////////////////////////////////////////////////////////////

vector<ModeEstimate*> Section2D::estimate_kz2_fourier()
{
  // Check values.

  if (    (global_section.leftwall  == no_wall)
       || (global_section.rightwall == no_wall) )
  {
    
    py_error("Function only defined for completely enclosed structures.");
    exit (-1);
  }

  // Calculate M and N.: TODO: better approximations

  // TODO: speed up by using more .reference.

  const Complex k0 = 2*pi/global.lambda;

  const Real W = floor(real(get_width())  * 1e10) / 1e10;
  const Real H = floor(real(get_height()) * 1e10) / 1e10;

  int M = int(sqrt(Real(M1*W/H)));
  int N = int(sqrt(Real(M1*H/W)));

  // Determine reflection coefficients of walls.

  Complex R_lower = slabs[0]->R_lower();
  Complex R_upper = slabs[0]->R_upper();

  Complex R_left  = global_section. leftwall == E_wall ? -1.0 : 1.0;
  Complex R_right = global_section.rightwall == E_wall ? -1.0 : 1.0;

  //if ( (abs(R_left+1.0) < 1e-6) && (abs(R_right+1.0) < 1e-6) ) // EE
  //    M--; // Note: don't do this for HH since it has a dummy solution.

  //if ( (abs(R_lower+1.0) < 1e-6) && (abs(R_upper+1.0) < 1e-6) ) // EE
  //    N--; // Note: don't do this for HH since it has a dummy solution.

  M = M1;
  N = M1;
  
  global_section.M = M;
  global_section.N = N;

  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;

  std::cout << "M N " << M << " " << N << " " << std::endl;

  Complex alpha0 = (abs( R_left*R_right-1.0) < 1e-10) ? 0 : pi/get_width()/2.;
  Complex  beta0 = (abs(R_lower*R_upper-1.0) < 1e-10) ? 0 : pi/get_height()/2.;

  if ((abs(R_lower-R_upper) > 1e-6) || ( abs(R_left-R_right) > 1e-6))
  {
    py_error("Error: only implemented for same wall type on opposite walls");
    exit(-1);
  }

  // Create eigenproblem.

  bool PML_present =  (    (abs(global_slab.lower_PML)    > 1e-12)
                        || (abs(global_slab.upper_PML)    > 1e-12)
                        || (abs(global_section.left_PML)  > 1e-12)
                        || (abs(global_section.right_PML) > 1e-12) );

  cMatrix F(2*MN,2*MN,fortranArray);
  cMatrix G(2*MN,2*MN,fortranArray);

  if (global_section.section_solver == NT)
    create_FG_NT(&F, & G, M, N, alpha0, beta0);
  else if (     (global_section.section_solver == L_anis) 
            || ((global_section.section_solver == L) && PML_present) )
    create_FG_li_biaxial(&F, & G, M, N, real(alpha0), real(beta0));
  else if (global_section.section_solver == L)
    create_FG_li(&F, & G, M, N, alpha0, beta0);
   else if (global_section.section_solver == ASR_2D 
         || global_section.section_solver == ASR_2D_stretched)
  	create_FG_ASR(&F, & G, M, N, alpha0, beta0);

  cMatrix FG(2*MN,2*MN,fortranArray);  
  FG.reference(multiply(F,G));

  bool reduced = global_section.reduced_eigenmatrix;
  
  if (!reduced)
  {
    cVector E(2*MN,fortranArray); 
    cMatrix eig(2*MN,2*MN,fortranArray); 

    if (global.stability == normal)
      E = eigenvalues(FG, &eig);
    else
      E = eigenvalues_x(FG, &eig);

    //for (int i=1; i<=eig.columns(); i++)
    //{
    //  std::cout << "eig " << i << " " << sqrt(E(i)/k0/k0)/k0 << std::endl;
    
    //for (int m=-M; m<=M; m++)
    // for (int n=-N; n<=N; n++)
    // {
    //    int i1 = (m+M+1) + (n+N)*(2*M+1);
    //    std::cout << i << " " <<m << " " << n << " " 
    //              << eig(i1,i)/eig(1,i) <<eig(MN+i1,i)/eig(1,i)<< std::endl;
    //  }
    //std::cout << std::endl;
    //}

    // The following code is only used to generate output.

    if (global_section.print_estimates)
    {
      // Sorting the E(i) vector might be a better idea.

      vector<Complex> E_temp;
      vector<Complex> E_final;

      for (int i=1; i <=E.rows(); i++)
        E_temp.push_back(E(i));

      std::sort(E_temp.begin(), E_temp.end(), RealSorter());

      vector<Complex> neff_temp;
      vector<Complex> neff;
      
      for (int i=1; i<=E_temp.size(); i++)
        neff_temp.push_back(sqrt(E_temp[i]/k0/k0)/k0);

      //std::sort(neff_temp.begin(), neff_temp.end(),RealSorter());

      Real n_eff_max = global_section.n_eff_max;
      int nov        = global_section.number_of_values;

      for (unsigned int i=0; i<neff_temp.size(); i++)
        if (real(neff_temp[i]) < n_eff_max)
        {
          neff.push_back(neff_temp[i]);
          E_final.push_back(E_temp[i]);
        }

      //std::sort(neff.begin() , neff.end(), RealSorter())

      for (int i=neff.size()-1; i>=neff.size()-1-nov; i--)
        std::cout << "n_eff_full " << neff.size()-i-1  
                  << " " << real(neff[i]) << " " << imag(neff[i]) 
                  << " " << "E " << real(E_final[i]) << " " << imag(E_final[i])
                  << std::endl;
    }
  }
   
  //
  // Reduced eigenvalue problem.
  //

  int M_ = int((m_+1)/2);  
  int N_ = int((n_+1)/2);

  int MN_ = M_*N_;

  cMatrix FG_(2*MN_,2*MN_,fortranArray);

  IndexConvertor f(M,N);

  vector<int> rows_to_keep;
  for (int m=0; m<=M; m++)
    for (int n=0; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);
      rows_to_keep.push_back(i1);      
      rows_to_keep.push_back(MN+i1);
    }
  
  std::sort(rows_to_keep.begin(), rows_to_keep.end());

  // Determine reduction coefficients.

  Real c[16];
  
  // Case 1: upper/lower=H, left/right=E

  if ((abs(R_lower-1.0) < 1e-6) && (abs(R_left+1.0) < 1e-6))
  {
    c[0]=1; c[1]=0;
  
    c[2]=1; c[3]=1; c[4]=0; c[5]=0;
    c[6]=1; c[7]=1; c[8]=0; c[9]=0;

    c[10]= 1; c[11]= 1; c[12]=1;
    c[13]=-1; c[14]=-1; c[15]=1;
  }
  
  // Case 2: upper/lower=E, left/right=H

  if ((abs(R_lower+1.0) < 1e-6) && (abs(R_left-1.0) < 1e-6))
  {
  c[0]=0; c[1]=1;
  
  c[2]=0; c[3]=0; c[4]=1; c[5]=1;
  c[6]=0; c[7]=0; c[8]=1; c[9]=1;

  c[10]=-1; c[11]=-1; c[12]=1;
  c[13]= 1; c[14]= 1; c[15]=1;
  }

  // Case 3: upper/lower=H, left/right=H

  if ((abs(R_lower-1.0) < 1e-6) && (abs(R_left-1.0) < 1e-6))
  {  
    c[0]=0; c[1]=0;
  
    c[2]=0; c[3]= 0; c[4]=1; c[5]=-1;
    c[6]=1; c[7]=-1; c[8]=0; c[9]=0;

    c[10]= 1; c[11]=-1; c[12]=-1;
    c[13]=-1; c[14]= 1; c[15]=-1;
  }
  
  // Case 4: upper/lower=E, left/right=E

  if ((abs(R_lower+1.0) < 1e-6) && (abs(R_left+1.0) < 1e-6))
  {  
    c[0]=0; c[1]=0;
  
    c[2]=1; c[3]=-1; c[4]=0; c[5]=0;
    c[6]=0; c[7]= 0; c[8]=1; c[9]=-1;

    c[10]=-1; c[11]= 1; c[12]=-1;
    c[13]= 1; c[14]=-1; c[15]=-1;
  }

  // Do actual reduction.
  
  for (int i1=1; i1<=2*MN_; i1++)
  {
    int i1_ = rows_to_keep[i1-1];

    // j == 0, l == 0

    FG_(i1,1)     = c[0]*FG(i1_,   f(0, 0));
    FG_(i1,MN_+1) = c[1]*FG(i1_,MN+f(0, 0));

    // j == 0, l != 0

    for (int l=1; l<=N; l++)
    {
      int i2 = (0+1) + l*(M+1);

      FG_(i1,i2)     = c[2]*FG(i1_,   f(0, l)) + c[3]*FG(i1_,   f(0, -l));
      FG_(i1,MN_+i2) = c[4]*FG(i1_,MN+f(0, l)) + c[5]*FG(i1_,MN+f(0, -l));
    }

    for (int j=1; j<=M; j++)
    {
     // j != 0, l == 0

      FG_(i1,j+1)     = c[6]*FG(i1_,   f(j, 0)) + c[7]*FG(i1_,   f(-j, 0));
      FG_(i1,MN_+j+1) = c[8]*FG(i1_,MN+f(j, 0)) + c[9]*FG(i1_,MN+f(-j, 0));

      // j != 0, l != 0

      for (int l=1; l<=N; l++)
      {
        int i2 = (j+1) + l*(M+1);

        FG_(i1,i2) =         FG(i1_,   f( j, l))
                     + c[10]*FG(i1_,   f( j,-l))
                     + c[11]*FG(i1_,   f(-j, l))
                     + c[12]*FG(i1_,   f(-j,-l));

        FG_(i1,i2+MN_) =     FG(i1_,MN+f( j, l))
                     + c[13]*FG(i1_,MN+f( j,-l))
                     + c[14]*FG(i1_,MN+f(-j, l))
                     + c[15]*FG(i1_,MN+f(-j,-l));
      }
    }
  }

  // Solve reduced eigenvalue problem.

  cVector E_(2*M_*N_,fortranArray);
  cMatrix eig_(2*M_*N_,2*M_*N_,fortranArray); 
  
  if (global.stability == normal)
    E_ = eigenvalues(FG_, &eig_);
  else
    E_ = eigenvalues_x(FG_, &eig_);

  // The following code is only used to generate output.

  if (global_section.print_estimates)
  {
    vector<Complex> E_temp;

    for (int i=1; i<=E_.rows(); i++)
      E_temp.push_back(E_(i));

    std::sort(E_temp.begin(), E_temp.end(), RealSorter());

    vector<Complex> neff_temp;

    for (int i=0; i<E_temp.size(); i++) // BUG PDB
      neff_temp.push_back(sqrt(E_temp[i]/k0/k0)/k0);

    int nov = global_section.number_of_values;

    vector<Complex> E_final;
    vector<Complex> neff;

    for (unsigned int i=0; i<neff_temp.size(); i++)
      if (real(neff_temp[i]) < global_section.n_eff_max)
      {
        neff.push_back(neff_temp[i]);
        E_final.push_back(E_temp[i]);
      }

    for (int i=neff.size()-1; i>=neff.size()-nov; i--)
      std::cout << "n_eff " << N << " " << neff.size()-i-1 << " " 
                << real(neff[i]) << " " << imag(neff[i]) << std::endl;
  }

  /*
  for (int i=1; i<=eig_.columns(); i++)
  {
    std::cout << "eig_ " << i << " " << sqrt(E_(i)/k0/k0)/k0 << std::endl;

    //for (int j=1; j<=eig_.rows(); j++)
    //  std::cout << j << " " << eig_(j,i) << std::endl;
    
    for (int m=0; m<=M; m++)
      for (int n=0; n<=N; n++)
      {
        int i1 = (m+1) + n*(M+1);
        std::cout << i << " " << " " <<m << " " << n << " " 
                  << eig_(i1,i) <<eig_(MN_+i1,i)<< std::endl;
      }
    std::cout << std::endl;
  }
  */
  
  vector<Complex> neff_;
  for (int i=1; i<= E_.rows(); i++)
    neff_.push_back(sqrt(E_(i)/k0/k0)/k0);
  std::sort(neff_.begin(), neff_.end(),RealSorter());

  //for (int i=0; i<neff_.size(); i++)
  //  std::cout << "reduced " << i << " " << neff_[i] << std::endl;

  std::cout << "Eigenproblem size " << 2*MN_ << std::endl;
  
  // Calculate field expansion.

  // TODO: try to do this using reduced matrices, or at least return
  // a reduced vector.
 
  cMatrix eig_big(2*MN,2*MN_,fortranArray); 

  for (int i1=1; i1<=2*MN_; i1++)
  {
    // j == 0, l == 0

    eig_big(   f(0, 0), i1) = c[0]*eig_(    1, i1);
    eig_big(MN+f(0, 0), i1) = c[1]*eig_(MN_+1, i1);

    // j == 0, l != 0

    for (int l=1; l<=N; l++)
    {
      int i2 = (0+1) + l*(M+1);

      eig_big(   f(0, l),i1) = c[2]*eig_(    i2,i1);
      eig_big(   f(0,-l),i1) = c[3]*eig_(    i2,i1);

      eig_big(MN+f(0, l),i1) = c[4]*eig_(MN_+i2,i1);
      eig_big(MN+f(0,-l),i1) = c[5]*eig_(MN_+i2,i1);
    }

    for (int j=1; j<=M; j++)
    {
     // j != 0, l == 0      

      eig_big(   f( j,0),i1) = c[6]*eig_(    j+1,i1);
      eig_big(   f(-j,0),i1) = c[7]*eig_(    j+1,i1);

      eig_big(MN+f( j,0),i1) = c[8]*eig_(MN_+j+1,i1);
      eig_big(MN+f(-j,0),i1) = c[9]*eig_(MN_+j+1,i1);

      // j != 0, l != 0

      for (int l=1; l<=N; l++)
      {
        int i2 = (j+1) + l*(M+1);

        eig_big(   f( j, l),i1) =       eig_(    i2,i1);
        eig_big(   f( j,-l),i1) = c[10]*eig_(    i2,i1);
        eig_big(   f(-j, l),i1) = c[11]*eig_(    i2,i1);
        eig_big(   f(-j,-l),i1) = c[12]*eig_(    i2,i1);

        eig_big(MN+f( j, l),i1) =       eig_(MN_+i2,i1);
        eig_big(MN+f( j,-l),i1) = c[13]*eig_(MN_+i2,i1);
        eig_big(MN+f(-j, l),i1) = c[14]*eig_(MN_+i2,i1);
        eig_big(MN+f(-j,-l),i1) = c[15]*eig_(MN_+i2,i1);
      }
    }
  }

/*
  for (int i=1; i<=eig_big.columns(); i++)
  {
    std::cout << "eig_big " << i << " " << sqrt(E_(i)/k0/k0)/k0 << std::endl;
    
    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
      {
        int i1 = (m+M+1) + (n+N)*(2*M+1);
        std::cout << i << " " <<m << " " << n << " " 
                  << eig_big(i1,i)/eig_big(1,i) 
                  <<eig_big(MN+i1,i)/eig_big(1,i)<< std::endl;
      }
    std::cout << std::endl;
  }
*/

  // Calculate H fields from E fields.

  cMatrix eig_big_H(2*MN,2*MN_,fortranArray);
  eig_big_H.reference(multiply(G,eig_big));

  // Return estimates.

  vector<ModeEstimate*> estimates;

  blitz::Range r1(1,MN); blitz::Range r2(MN+1,2*MN);

  bool TEM = false;
  
  //for (int i=1; i<=E_.rows(); i++)
  //  std::cout << E_(i)/k0/k0 << " " << sqrt(E_(i)/k0/k0)/k0 << std::endl;

  const Real Y0 = sqrt(eps0/mu0);
 
  for (int i=1; i<=E_.rows(); i++)
  { 
    // Compensate for spurious TEM mode.

    if (abs(E_(i) - pow(k0*k0*core->n(), 2)) < 1e-6)
    {
      TEM = true;
      continue;
    }
    
    if ((abs(E_(i)) > 1e-6))
    {
      Complex kz2 = E_(i)/k0/k0;
      Complex kz = sqrt(kz2);

      if (imag(kz) > 0)
        kz = -kz;

      if (abs(imag(kz)) < abs(real(kz)))
        if (real(kz) < 0)
          kz = -kz;

      // Compensate for numerical instability in the presence of PML.

      if ((real(kz) > imag(kz)) && (imag(kz) > 1e-8) 
          && (global_section.mode_correction == snap))
      {
        std::cout << "Snapping " << kz/2./pi*global.lambda
                  << " to the real axis" << std::endl;
        kz = real(kz);
      }
      
      cVector* Ex = new cVector(MN,fortranArray);
      cVector* Ey = new cVector(MN,fortranArray);
      cVector* Hx = new cVector(MN,fortranArray);
      cVector* Hy = new cVector(MN,fortranArray);

      *Ex = eig_big  (r1,i); 
      *Ey = eig_big  (r2,i);
      *Hx = eig_big_H(r1,i)/kz/k0*Y0;
      *Hy = eig_big_H(r2,i)/kz/k0*Y0;

      Ex->makeUnique(); Ey->makeUnique();
      Hx->makeUnique(); Hy->makeUnique();

      ModeEstimate* est = new ModeEstimate(kz*kz, Ex,Ey, Hx,Hy);
      estimates.push_back(est);
    }  
  }

  std::sort(estimates.begin(), estimates.end(), kz2_sorter());

  //for (int i=0; i<estimates.size(); i++)
  //  std::cout << "sorted " << i << " " << sqrt(estimates[i]->kz2)/k0 
  //            << std::endl;

  if (estimates.size() > 4)
    for (int i=0; i<(TEM ? 3 : 4); i++)
    {
      delete estimates.back();
      estimates.erase(estimates.end()-1);
    }

  std::ostringstream s;
  s << "Found " << estimates.size() << " estimates.";
  py_print(s.str());

  return estimates;
}
  


/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_from_estimates()
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes_from_estimates()
{
  // Check M1.

  int n = M1/2;
  
  if (2*n != M1)
    py_print("Warning: changing M1 to even number.");

  M1 = 2*n;

  // Find min, max and cladding eps mu.

  vector<Complex> eps_mu;
  for (unsigned int i=0; i<materials.size(); i++)
    eps_mu.push_back(materials[i]->eps_mu());

  remove_copies(&eps_mu, 1e-23);
  std::sort(eps_mu.begin(), eps_mu.end(), RealSorter());
  
  Complex min_eps_mu, max_eps_mu, clad_eps_mu;
  
  if (eps_mu.size() >= 2)
  {
     min_eps_mu = eps_mu[0];
     max_eps_mu = eps_mu[eps_mu.size()-1];
    clad_eps_mu = eps_mu[eps_mu.size()-2];
  }
  else
    min_eps_mu = max_eps_mu = clad_eps_mu = eps_mu[0];

  const Complex C0 = pow(2*pi/global.lambda, 2) / eps0 / mu0;

  // Get estimates.
                                  
  vector<ModeEstimate*> estimates;
                                   
  if (user_estimates.size() != 0) // User provided estimates.
  {
    for (unsigned int i=0; i<user_estimates.size(); i++)
    {
      ModeEstimate* est 
        = new ModeEstimate(pow(2*pi/global.lambda*user_estimates[i], 2));
      estimates.push_back(est);
    } 
  }
  else
  {
    Complex max_eps_eff;
   
    if (real(max_eps_mu)*real(min_eps_mu) < 0.0) // Surface plasmon.
        max_eps_eff = 1.0/(1.0/max_eps_mu + 1.0/min_eps_mu);
    else
      max_eps_eff = max_eps_mu;

    if (abs(real(max_eps_eff)) < abs(real(max_eps_mu)))
      max_eps_eff = max_eps_mu;

    Real max_kz = abs(2*pi/global.lambda*sqrt(max_eps_eff/eps0/mu0));
    vector<ModeEstimate*> estimates_0;   

    if (global_section.section_solver == OS)
      estimates_0 = estimate_kz2_omar_schuenemann();  
    else
      estimates_0 = estimate_kz2_fourier();

    Real C = global_section.estimate_cutoff;
    for (unsigned int i=0; i<estimates_0.size(); i++)
      if (    (real(sqrt(estimates_0[i]->kz2)) < C*max_kz) 
           || (global_section.keep_all_estimates == true))
          estimates.push_back(estimates_0[i]);
  }

  if (global_section.keep_all_estimates == true)
  {    
    global.N = estimates.size();
    std::ostringstream s;
    s << "Setting N to " << estimates.size() << ".";
    py_print(s.str());
  }

  while (estimates.size() > global.N)
  {
    delete estimates.back();
    estimates.erase(estimates.end()-1);
  }

  // Split off modes to be corrected.

  //bool metallic = (real(min_eps_mu) < 0);
  bool metallic = global_slab.low_index_core;

  Complex threshold_eps_mu = metallic ? max_eps_mu : clad_eps_mu;
  Complex        kt_eps_mu = metallic ? max_eps_mu :  min_eps_mu;

  vector<Complex> kt2_coarse;

  for (unsigned int i=0; i<estimates.size(); i++)
  {
    Complex kz = sqrt(estimates[i]->kz2);

    // Mode to be corrected.

    if (    (user_estimates.size() != 0)
         || ( global_section.mode_correction == full)
         || (    (global_section.mode_correction == guided_only)
              && ( real(kz/2./pi*global.lambda) >
                   real(sqrt(threshold_eps_mu/eps0/mu0)))))
    {
      Complex kt2 = C0*kt_eps_mu - estimates[i]->kz2;
      kt2_coarse.push_back(kt2);
    }
    else // Uncorrected mode.
    {
      if (imag(kz) > 0)
        kz = -kz;

      if (abs(imag(kz)) < abs(real(kz)))
        if (real(kz) < 0)
          kz = -kz;

      Section2D_Mode* newmode = new Section2D_Mode(TE_TM, kz, this,
                                estimates[i]->Ex, estimates[i]->Ey,
                                estimates[i]->Hx, estimates[i]->Hy, 
                                false);
      newmode->normalise();
      modeset.push_back(newmode);
    }
  }
  
  user_estimates.clear();

  // Refine modes using transcendental function.

  if (kt2_coarse.size())
  {
    if (global_section.mode_correction == guided_only)
    {
      std::ostringstream s;
      s << "Refining estimates with n_eff higher than " 
        << sqrt(clad_eps_mu/eps0/mu0) << "...";
      py_print(s.str());
    }
    else
      py_print("Refining estimates...");
  }

  kt2_to_neff transform(C0*kt_eps_mu);
  
  SectionDisp disp(left, right, global.lambda, M2, symmetric);

  //for (Real n_eff=2.2; n_eff<=2.25; n_eff+=0.001)
  //{ 
  //    for (Real n_eff_imag=0.; n_eff_imag<=0.1; n_eff_imag +=0.001)
  //    {
  //      Complex kt = sqrt(C0*min_eps_mu-pow(2.*pi/global.lambda
  //                       *(n_eff-n_eff_imag*I),2));
  //      //std::cout << n_eff << " " << real(kt) << " " << imag(kt) 
  //      //      << " " << abs(disp(kt)) << std::endl;
  //      std::cout << n_eff << " " << n_eff_imag << " " 
  //            << " " << abs(disp(kt)) << std::endl;
  //  }
  //}
  // exit(-1);

  vector<Complex> kt2 = mueller(disp, kt2_coarse, 1e-8, 100, &transform, 2);

  // Eliminate false zeros.

  for (unsigned int k=0; k<slabs.size(); k++)
    for (unsigned int i=1; i<=slabs[k]->N(); i++)
    {
      Complex beta_i = slabs[k]->get_mode(i)->get_kz();
      Complex kt2_i = C0*kt_eps_mu - beta_i*beta_i;

      remove_elems(&kt2,  kt2_i, 1e-6);
      remove_elems(&kt2, -kt2_i, 1e-6);
    }

  // Add corrected modes to modeset.

  py_print("Creating mode profiles...");

  for (unsigned int i=0; i<kt2.size(); i++)
  {
    Complex kz = sqrt(C0*kt_eps_mu - kt2[i]);
    
    if (real(kz) < 0) 
      kz = -kz;

    if (abs(real(kz)) < 1e-12)
      if (imag(kz) > 0)
        kz = -kz;

    Section2D_Mode* newmode = new Section2D_Mode(TE_TM, kz, this);

    newmode->normalise();

    modeset.push_back(newmode);
  }

  // Gracefully handle case when not enough modes were found.

  if (modeset.size() < global.N)
  {
    std::ostringstream s;
    s << "Error: didn't find enough modes ("
      << modeset.size() << "/" << global.N << "). ";
    py_error(s.str());
    
    while (modeset.size() != global.N)
    {    
      Section2D_Mode* dummy = new Section2D_Mode(TE_TM, 0.0, this);
      modeset.push_back(dummy);
    }
  }

  // Free estimates.

  for (int i=0; i<estimates.size(); i++)
    delete estimates[i];

  sort_modes();

  py_print("Done.");

  return;

  // Test orthogonality.

  cMatrix O12(modeset.size(), modeset.size(), fortranArray);
  cMatrix O21(modeset.size(), modeset.size(), fortranArray);
  cMatrix O11(modeset.size(), modeset.size(), fortranArray);
  cMatrix O22(modeset.size(), modeset.size(), fortranArray);
  calc_overlap_matrices(this, &O12, &O21, &O11, &O22);
  std::cout << O11 << std::endl;
  std::cout << O12 << std::endl;

  //return;

  for (unsigned int i=0; i<modeset.size(); i++)
    for (unsigned int j=0; j<modeset.size(); j++)
        std::cout << i << " " << j << " " 
              << overlap_numeric(dynamic_cast<Section2D_Mode*>(modeset[i]),
                         dynamic_cast<Section2D_Mode*>(modeset[j]))
                  << std::endl << std::flush;

  exit(-1);
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
    
    Section2D_Mode *newmode = new Section2D_Mode(TE_TM,kz,this);
      
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

  if (    (global_section.leftwall  == no_wall)
       || (global_section.rightwall == no_wall) )
  {
    
    py_error("Function only defined for completely enclosed structures.");
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
          = new Section1D_Mode(TE_TM, kz, mode, this);

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
    last_gain_mat_n = global.gain_mat->n();
}

