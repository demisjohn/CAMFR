
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
  if ( (M1 < M2) || (M2 < global.N) )
    py_print("Warning: M1 > M2 > N not fulfilled.");

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
  right_ex.add_term(Term(*rightwall_sc));

  // Create left hand side expression.

  Expression left_ex;

  Slab* slab0 = dynamic_cast<Slab*>(ex.get_term(max_eps_i)->get_wg());
  Complex d0 = ex.get_term(max_eps_i)->get_d()/2.;
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
  left_ex.add_term(Term(*leftwall_sc));

  // Create Section.

  s = new Section2D(left_ex, right_ex, M1, M2);

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
    find_modes_from_series();

  // Remember wavelength and gain these modes were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}



/////////////////////////////////////////////////////////////////////////////
//
// Small function objects
//
/////////////////////////////////////////////////////////////////////////////
  
struct index_sorter
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
    {
      return ( real(beta_a) > real(beta_b) );
    }
};

struct loss_sorter
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
    {
      return ( abs(imag(sqrt(beta_a))) < abs(imag(sqrt(beta_b))) );
    }
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
// Section2D::estimate_kz2
//
/////////////////////////////////////////////////////////////////////////////

cVector Section2D::estimate_kz2()
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
  if (global.stability == normal)
    kz2 = eigenvalues(E);
  else
    kz2 = eigenvalues_x(E);

  return kz2;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier
//
//  Returns fourier expansion of order m = -M,...,M of a piecewise 
//  continuous function f given by:
//
//    disc[0] -> disc[1]   : f[0]
//    ...
//    disc[n] -> disc[n+1] : f[n]
//
//  Basis functions are exp(j.m.2.pi/d.x). (d can be overridden).
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier(const vector<Complex>& f, const vector<Complex>& disc, int M,
                const Complex* d=0)
{
  const Complex D = d ? *d : disc.back() - disc.front();
  const Complex K = 2.*pi/D;

  cVector result(2*M+1,fortranArray);

  for (int m=-M; m<=M; m++)
  {
    Complex result_m = 0.0;

    for (unsigned int k=0; k<int(disc.size()-1); k++)
    {
      Complex factor;
      if (m==0)
        factor = disc[k+1]-disc[k];
      else
      {
        const Complex t = -I*Real(m)*K;
        factor = (exp(t*disc[k+1])-exp(t*disc[k])) / t;
      }

      result_m += factor * f[k];
    }
    
    result(m+M+1) = result_m / D;
  }

  if (abs(D) < 1e-12)
    result = 0.0;
  
  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
//  Returns 2D fourier expansion of order m_x = -M,...,M and n_y = -N,...,N 
//  of dielectric profile give by a sequence of slabs:
//
//    disc[0] -> disc[1]   : slabs[0]
//    ...
//    disc[n] -> disc[n+1] : slabs[n]
//
//  Basis functions are exp(j.m.2.pi/d.x).exp(j.n.2.pi/d.y).
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_eps_2D(const vector<Slab*>& slabs, 
                       const vector<Complex>& disc, int M, int N, 
                       bool invert=false)
{
  const Complex Lx = disc.back() - disc.front();
  cMatrix result(2*M+1,2*N+1,fortranArray);
  result = 0.0;
  
  for (int i=0; i<slabs.size(); i++)
  {
    // Calculate 1D fourier transform of eps(y) profile.

    vector<Complex> disc_i_y(slabs[i]->get_discontinuities());
    disc_i_y.insert(disc_i_y.begin(), 0.0);
    
    vector<Complex> f_i_y;
    for (int k=0; k<disc_i_y.size()-1; k++)
    {
      Complex eps = slabs[i]->eps_at(Coord(disc_i_y[k],0,0,Plus))/eps0;
      
      if (invert)
        eps = 1.0/eps;
      
      f_i_y.push_back(eps);
    }
    
    cVector fourier_1D_y(2*N+1,fortranArray);   
    fourier_1D_y = fourier(f_i_y, disc_i_y, N);

    // Calculate pseudo 1D fourier transform in x direction.

    vector<Complex> disc_i_x;
    disc_i_x.push_back(disc[i]);
    disc_i_x.push_back(disc[i+1]);

    vector<Complex> f_i_x; 
    f_i_x.push_back(1.0);

    cVector fourier_1D_x(2*M+1,fortranArray);
    fourier_1D_x = fourier(f_i_x, disc_i_x, M, &Lx);

    // Calculate 2D fourier transform.

    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
        result(m+M+1,n+N+1) += fourier_1D_x(m+M+1) * fourier_1D_y(n+N+1);
  }
  
  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_from_series()
//
/////////////////////////////////////////////////////////////////////////////

cVector Section2D::estimate_kz2_li()
{
#if 0
  vector<Complex> disc;  
  vector<Complex> f;

  disc.push_back(0.0);
  f.push_back(0.0);
  disc.push_back(1.0);
  f.push_back(1.0);  
  disc.push_back(2.0);
  f.push_back(0.0);  
  disc.push_back(4.0);
  
  int M = 10;
  cVector F(2*M+1, fortranArray);
  F = fourier(f,disc,M);

  for (Real x=0; x<=4; x+=.1)
  {
    Complex result = 0.0;
    for (int m=-M; m<=M; m++)
      result += F(m+M+1)*exp(I*Real(m)*2.*pi/4.*x);

    std::cout << x << " " << real(result) 
              << " " << imag(result) << std::endl;
  }
#endif

#if 0
  vector<Complex> disc(discontinuities);
  disc.insert(disc.begin(), 0.0);

  int M = 20;
  int N = 20;
  
  cMatrix F(2*M+1,2*N+1,fortranArray);
  F = fourier_eps_2D(slabs, disc, M, N);

  Complex d_x = disc.back();
  Complex d_y = slabs[0]->get_width();

  for (Real x=0; x<=real(d_x); x+=.05)
  {
    for (Real y=0; y<=real(d_y); y+=.05)
    {
      Complex result = 0.0;
      for (int m=-M; m<=M; m++)      
        for (int n=-N; n<=N; n++)
          result += F(m+M+1,n+N+1)*exp(I*Real(m)*2.*pi/d_x*x)
                                  *exp(I*Real(n)*2.*pi/d_y*y);

      std::cout << x << " " << y << " " << real(result) << std::endl;
    }
  }
#endif

  // Calculate M and N. 
  // TODO: take aspect ratio into account.

  const Complex k0 = 2*pi/global.lambda;

  int M = int(sqrt(Real(M1)));
  int N = int(sqrt(Real(M1)));

  // Calculate eps and inv_eps matrix.

  vector<Complex> disc(discontinuities);
  disc.insert(disc.begin(), 0.0);
  
  cMatrix eps_(4*M+1,4*N+1,fortranArray);
  eps_ = fourier_eps_2D(slabs, disc, 2*M, 2*N);

  cMatrix inv_eps_(4*M+1,4*N+1,fortranArray);
  inv_eps_ = fourier_eps_2D(slabs, disc, 2*M, 2*N, true);  

  const int MN = (2*M+1)*(2*N+1);
  
  cMatrix     eps(MN,MN,fortranArray);
  cMatrix inv_eps(MN,MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {
      int i1 = (m+M+1) + (n+N)*(2*M+1);
      
      for (int j=-M; j<=M; j++)
        for (int l=-N; l<=N; l++)
        {
          int i2 = (j+M+1) + (l+N)*(2*M+1);

              eps(i1,i2) =     eps_(m-j + 2*M+1, n-l + 2*N+1);          
          inv_eps(i1,i2) = inv_eps_(m-j + 2*M+1, n-l + 2*N+1);
        }
    }
  
  // Calculate alpha and beta vector.

  Complex alpha0, beta0;
  alpha0 = beta0 = 0.0;

  cVector alpha(MN,fortranArray);
  cVector  beta(MN,fortranArray);

  for (int m=-M; m<=M; m++)
    for (int n=-N; n<=N; n++)
    {      
      int i = (m+M+1) + (n+N)*(2*M+1);
      alpha(i) = alpha0 + m*2.*pi/get_width();
       beta(i) =  beta0 + n*2.*pi/get_height();     
    }

  // Constuct F matrix.

  cMatrix F(2*MN,2*MN,fortranArray);
  
  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      F(i1,   i2)    =  alpha(i1) * inv_eps(i1,i2) *  beta(i2);      
      F(i1,   i2+MN) = -alpha(i1) * inv_eps(i1,i2) * alpha(i2);        
      F(i1+MN,i2)    =   beta(i1) * inv_eps(i1,i2) *  beta(i2);
      F(i1+MN,i2+MN) =  -beta(i1) * inv_eps(i1,i2) * alpha(i2);      

      if (i1==i2)   
      {
        F(i1,   i2+MN) += k0*k0;
        F(i1+MN,i2)    -= k0*k0;
      }
    }
  
  // Construct G matrix.

  cMatrix G(2*MN,2*MN,fortranArray);
  
  for (int i1=1; i1<=MN; i1++)
    for (int i2=1; i2<=MN; i2++)
    {
      G(i1,   i2)    = (i1==i2) ? -alpha(i1)*beta(i2) : 0.0;
      G(i1,   i2+MN) = -k0*k0 * eps(i1,i2);        
      G(i1+MN,i2)    =  k0*k0 * eps(i1,i2);
      G(i1+MN,i2+MN) = (i1==i2) ?  alpha(i1)*beta(i2) : 0.0;

      if (i1==i2) 
      {     
        G(i1,   i2+MN) += alpha(i1)*alpha(i1);  
        G(i1+MN,i2)    -=  beta(i2)*beta(i2);
      }
      
    }
  
  // Solve eigenproblem.

  cMatrix FG(2*MN,2*MN,fortranArray); 
  FG.reference(multiply(F,G));

  cVector E(2*MN,fortranArray);
  
  if (global.stability == normal)
    E = eigenvalues(FG);
  else
    E = eigenvalues_x(FG);

  E /= k0*k0;

  std::cout << "Eigenproblem size " << 2*MN << std::endl;
  for (int i=1; i<2*MN; i++)
    std::cout << i << " " << E(i) << " " << sqrt(E(i)) << " " 
              << sqrt(E(i))/2./pi*global.lambda << std::endl;
  
  return E;
}



/////////////////////////////////////////////////////////////////////////////
//
// Section2D::find_modes_from_series()
//
/////////////////////////////////////////////////////////////////////////////

void Section2D::find_modes_from_series()
{
  // Check M1.

  int n = M1/2;
  
  if (2*n != M1)
    py_print("Warning: changing M1 to even number.");

  M1 = 2*n;

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

  const Real C0 = pow(2*pi/global.lambda, 2) / eps0 / mu0;

  // Get estimates.

  vector<Complex> kz2_coarse;

  if (user_estimates.size() != 0) // User provided estimates
  {
    for (unsigned int i=0; i<user_estimates.size(); i++)
      kz2_coarse.push_back(pow(2*pi/global.lambda*user_estimates[i], 2));
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

    Real max_kz = real(2*pi/global.lambda*sqrt(max_eps_eff/eps0/mu0));
    
    cVector kz2(M1,fortranArray);

    if (global.eigen_calc == arnoldi)
      kz2 = estimate_kz2_li();
    else
      kz2 = estimate_kz2();

    //for (unsigned int i=1; i<=kz2.size(); i++)
    //  std::cout << "raw " << i << " "
    //            << sqrt(kz2(i))/2./pi*global.lambda << std::endl;

    for (unsigned int i=1; i<=M1; i++)
      if (real(sqrt(kz2(i))) < 1.01*max_kz)
        kz2_coarse.push_back(kz2(i));
  }
  user_estimates.clear();

  // Refine estimates.

  py_print("Refining estimates...");

  if (sort == highest_index)
    std::sort(kz2_coarse.begin(), kz2_coarse.end(), index_sorter());
  else
    std::sort(kz2_coarse.begin(), kz2_coarse.end(), loss_sorter());

  //for (unsigned int i=0; i<kz2_coarse.size(); i++)
  //  std::cout << i << " " << sqrt(kz2_coarse[i])/2./pi*global.lambda 
  //            << std::endl;  

  if (kz2_coarse.size() > global.N+5)
    kz2_coarse.erase(kz2_coarse.begin()+global.N+5, kz2_coarse.end());

  vector<Complex> kt_coarse;
  for (unsigned int i=0; i<kz2_coarse.size(); i++)
  {
    Complex kt = sqrt(C0*min_eps_mu - kz2_coarse[i]);

    if (imag(kt) < 0) 
      kt = -kt;

    if (abs(imag(kt)) < 1e-12)
      if (real(kt) > 0)
        kt = -kt;

    if ((abs(real(kt)) < .001) && (real(kt) < 0))
      kt -= 2*real(kt);

    kt_coarse.push_back(kt);
  }

  kt_to_neff transform(C0*min_eps_mu);
  SectionDisp disp(left, right, global.lambda, M2, symmetric);
  vector<Complex> kt = mueller(disp, kt_coarse, 1e-8, 100, &transform, 2);

  f = new SectionDisp(left, right, global.lambda, M2, symmetric); // TMP

  // Eliminate false zeros.

  for (unsigned int k=0; k<slabs.size(); k++)
    for (unsigned int i=1; i<=slabs[k]->N(); i++)
    {
      Complex beta_i = slabs[k]->get_mode(i)->get_kz();
      Complex kt_i = sqrt(C0*min_eps_mu - beta_i*beta_i);

      remove_elems(&kt,  kt_i, 1e-6);
      remove_elems(&kt, -kt_i, 1e-6);
    }

  // Create modeset.

  py_print("Creating mode profiles...");

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

  py_print("Done.");

  // Test orthogonality.

  return;

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
