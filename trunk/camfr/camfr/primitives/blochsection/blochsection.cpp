
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsection.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20050517
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include "blochsection.h"
#include "blochsectionmode.h"
#include "blochsectionoverlap.h"
#include "../slab/isoslab/slab.h"
#include "../section/section.h"
#include "../../math/calculus/fourier/fourier.h"

using std::vector;
using std::cout;
using std::endl;

#include "../../util/vectorutil.h"
#include "../../util/index.h"

/////////////////////////////////////////////////////////////////////////////
//
// BlochSectionGlobal
//
/////////////////////////////////////////////////////////////////////////////

BlochSectionGlobal global_blochsection = {0,0,0.0,0.0,0.0,0.0,0.5};



/////////////////////////////////////////////////////////////////////////////
//
// BlochSectionImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

void BlochSectionImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I, cMatrix* O_II_II)
{
  BlochSectionImpl* medium_I  = dynamic_cast<BlochSectionImpl*>(this);
  BlochSectionImpl* medium_II = dynamic_cast<BlochSectionImpl*>(w);

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
 
  // Overlap.

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
        (dynamic_cast<BlochSectionMode*>(medium_I ->get_mode(i)),
         dynamic_cast<BlochSectionMode*>(medium_II->get_mode(j)));

      (*O_II_I)(i,j) = overlap
        (dynamic_cast<BlochSectionMode*>(medium_II->get_mode(i)),
         dynamic_cast<BlochSectionMode*>(medium_I ->get_mode(j)));

      if (O_I_I) (*O_I_I)(i,j) = overlap
        (dynamic_cast<BlochSectionMode*>(medium_I ->get_mode(i)),
         dynamic_cast<BlochSectionMode*>(medium_I ->get_mode(j)));

      if (O_II_II) (*O_II_II)(i,j) = overlap
        (dynamic_cast<BlochSectionMode*>(medium_II->get_mode(i)),
         dynamic_cast<BlochSectionMode*>(medium_II->get_mode(j)));
    }
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection::BlochSection
//
/////////////////////////////////////////////////////////////////////////////

BlochSection::BlochSection(const Term& t)
{ 
  Slab* slab = dynamic_cast<Slab*>(t.get_wg());

  if ( !slab || !slab->is_uniform() )
  {
    py_error("Error: expected a uniform slab to initialise a section.");
    return;
  } 

  Complex W = t.get_d()
    + I*global_section.left_PML + I*global_section.right_PML;
  Complex H = slab->get_width() 
    + I*global_slab.lower_PML + I*global_slab.upper_PML;

  s = new UniformBlochSection(slab->get_core(), W, H);

  uniform = true;
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection::BlochSection
//
/////////////////////////////////////////////////////////////////////////////

BlochSection::BlochSection(Expression& expression) 
{
  // Uniform Section
  
  if (expression.get_size() == 1)
  { 
    Slab* slab = dynamic_cast<Slab*>(expression.get_term(0)->get_wg());

    if (!slab || !slab->is_uniform())
    {
      py_error("Error: expected a uniform slab to initialise a section.");
      return;
    } 

    Complex W = expression.get_term(0)->get_d()
      + I*global_section.left_PML + I*global_section.right_PML;
    Complex H = slab->get_width() 
      + I*global_slab.lower_PML + I*global_slab.upper_PML;

    s = new UniformBlochSection(slab->get_core(), W, H);

    uniform = true;
    core = s->get_core();

    return;
  }

  // Non-uniform BlochSection.

  s = new BlochSection2D(expression);
  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::BlochSection2D
//
/////////////////////////////////////////////////////////////////////////////

BlochSection2D::BlochSection2D(Expression& ex) : st(ex)
{
  // Determine core.

  uniform = false;

  materials = st.get_materials();

  unsigned int core_index = 0;
  for (unsigned int i=0; i<materials.size(); i++)
    if ( real(materials[i]->eps_mu()) > real(materials[core_index]->eps_mu()) )
      core_index = i;

  core = materials[core_index];

  // Determine slabs.

  Complex z = 0;

  Expression ex_flat = ex.flatten();
  for (int i=ex_flat.get_size()-1; i>=0; i--)
  {
    const Term* t = ex_flat.get_term(i);
    if (t->get_type() != WAVEGUIDE)
      continue;
    
    z += t->get_d();

    discontinuities.push_back(z);
    slabs.push_back(dynamic_cast<Slab*>(t->get_wg()));
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::BlochSection2D
//
/////////////////////////////////////////////////////////////////////////////

BlochSection2D::BlochSection2D(const BlochSection2D& section)
{
  st        = section.st;
  core      = section.core;
  uniform   = section.uniform;
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::operator=
//
/////////////////////////////////////////////////////////////////////////////

BlochSection2D& BlochSection2D::operator=(const BlochSection2D& section)
{
  if (this == &section)
    return *this;

  st        = section.st;
  core      = section.core;
  uniform   = section.uniform;
  
  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::no_gain_present
//
/////////////////////////////////////////////////////////////////////////////

bool BlochSection2D::no_gain_present() const
{
  return st.no_gain_present();
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::material_at
//
/////////////////////////////////////////////////////////////////////////////

Material* BlochSection2D::material_at(const Coord& c) const
{
  Complex c1 = c.c1;
  Limit c1_limit = c.c1_limit;

  c1 = st.get_total_thickness() - c.c1;
  c1_limit = (c.c1_limit == Plus) ? Min : Plus;
  
  Coord c_new(c.c2, c.z, c1, c.c2_limit, c.z_limit, c1_limit);

  return st.material_at(c_new);
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



/////////////////////////////////////////////////////////////////////////////
//
// rotate_slabs_()
//
//  Given a sequence of slabs in the x direction, transform it to
//  a sequence of slabs in the y direction.
//
/////////////////////////////////////////////////////////////////////////////

// TODO: remove redundant copy

void rotate_slabs_(const vector<Slab*>& slabs, const vector<Complex>& disc,
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
// BlochSection2D::create_FG_NT
//
//   Noponen/Turunen formulation
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D::create_FG_NT(cMatrix* F, cMatrix* G, int M, int N,
                                  const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2*pi/global.lambda;

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
  eps_ = fourier_2D(disc_x, disc_y, f_eps, 2*M, 2*N);

  cMatrix inv_eps_(4*M+1,4*N+1,fortranArray);
  inv_eps_ = fourier_2D(disc_x, disc_y, f_inv_eps, 2*M, 2*N);
  
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

      alpha(i) = alpha0 + m*2.*pi/get_width();
       beta(i) =  beta0 + n*2.*pi/get_height();
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
// BlochSection2D::create_FG_li
//
//  Li formulation
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D::create_FG_li(cMatrix* F, cMatrix* G, int M, int N,
                                  const Complex& alpha0, const Complex& beta0)
{
  const Complex k0 = 2.0*pi/global.lambda;

  // Construct data structures for fourier analysis.

  vector<Complex> disc_x(discontinuities);
  disc_x.insert(disc_x.begin(), 0.0);

  vector<vector<Complex> > disc_y, f_eps;

  for (int i=0; i<disc_x.size()-1; i++)
  {
    vector<Complex> disc_y_i(slabs[i]->get_discontinuities());
    disc_y_i.insert(disc_y_i.begin(), 0.0);

    vector<Complex> f_eps_i;
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
  rotate_slabs_(slabs, disc_x, &slabs_rot, &disc_x_rot);

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

  eps_ = fourier_2D(disc_x, disc_y, f_eps, 2*M, 2*N);
  
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
  t.reference(fourier_2D_split(disc_x_rot,disc_y_rot,f_eps_rot,N,M));
  if (global.stability != SVD)
    eps_y_x.reference(invert(t));
  else
    eps_y_x.reference(invert_svd(t));

  eps_x_y.reference(fourier_2D_split(disc_x_rot,disc_y_rot,f_inv_eps_rot,N,M));
  
  // Calculate alpha and beta vector.

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

  // Free rotated slabs.

  for (int i=0; i<slabs_rot.size(); i++)
    delete slabs_rot[i];
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::create_FG_li_biaxial
//
//  Li formulation for biaxial media
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D::create_FG_li_biaxial(cMatrix* F, cMatrix* G, 
                                          int M, int N,
                                          const Complex& alpha0, 
                                          const Complex& beta0)
{
  const Complex k0 = 2.*pi/global.lambda;

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

  eps_1. reference(fourier_2D_split(disc_x, disc_y, f_eps_1, M, N));
   mu_1. reference(fourier_2D_split(disc_x, disc_y,  f_mu_1, M, N)); 
 
  eps_2_.reference(fourier_2D_split(disc_x, disc_y, f_eps_2, M, N));
   mu_2_.reference(fourier_2D_split(disc_x, disc_y,  f_mu_2, M, N));  

  cMatrix eps_3_f(4*M+1,4*N+1,fortranArray); 
  cMatrix  mu_3_f(4*M+1,4*N+1,fortranArray);
  
  eps_3_f = fourier_2D(disc_x, disc_y, f_eps_3, 2*M, 2*N);
   mu_3_f = fourier_2D(disc_x, disc_y,  f_mu_3, 2*M, 2*N);

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

      alpha(i) = alpha0 + m*2.*pi/real(get_width());
       beta(i) =  beta0 + n*2.*pi/real(get_height());
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
// BlochSection2D::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D::find_modes()
{  
  // Check values.

  if (real(global.lambda) == 0)
  {
    py_error("Error: wavelength not set.");
    return;
  }
  
  // TODO: improve check for recalc needed.
  
  //if (!recalc_needed())
  //  return;

  // TODO: speed up by using more .reference.

  const Complex k0 = 2*pi/global.lambda;
  
  int M = global_blochsection.Mx;
  int N = global_blochsection.My;

  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;

  Complex alpha0 = global_blochsection.alpha0;
  Complex  beta0 = global_blochsection.beta0;

  // Create eigenproblem.

  bool PML_present =  (    (abs(global_slab.lower_PML)         > 1e-12)
                        || (abs(global_slab.upper_PML)         > 1e-12)
                        || (abs(global_blochsection.left_PML)  > 1e-12)
                        || (abs(global_blochsection.right_PML) > 1e-12) );

  cMatrix F(2*MN,2*MN,fortranArray);
  cMatrix G(2*MN,2*MN,fortranArray);

  if ( (global_section.section_solver == L_anis) || (PML_present) )
    create_FG_li_biaxial(&F, & G, M, N, real(alpha0), real(beta0));
  else
    create_FG_li(&F, & G, M, N, alpha0, beta0);

  //create_FG_NT(&F, & G, M, N, alpha0, beta0);

  cMatrix FG(2*MN,2*MN,fortranArray);  
  FG.reference(multiply(F,G));

  cVector E(2*MN,fortranArray); 
  cMatrix eig(2*MN,2*MN,fortranArray); 

  if (global.stability == normal)
    E = eigenvalues(FG, &eig);
  else
    E = eigenvalues_x(FG, &eig);

  //vector<Complex> neff;
  //for (int i=1; i<= E.rows(); i++)
  //  neff.push_back(sqrt(E(i)/k0/k0)/k0);
  //std::sort(neff.begin(), neff.end(), RealSorter());

  //std::cout << "M N " << M << " " << N << std::endl;
  //std::cout << "Eigenproblem size " << 2*MN << std::endl;
  //for (int i=neff.size()-1; i>=0; i--)
  //  std::cout << "n_eff " << neff.size()-i-1  << " " << neff[i] << std::endl;

  // Calculate H fields from E fields.

  cMatrix eig_H(2*MN,2*MN,fortranArray);
  eig_H.reference(multiply(G,eig));

  // Clear modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Construct modes.

  blitz::Range r1(1,MN); blitz::Range r2(MN+1,2*MN);

  const Real Y0 = sqrt(eps0/mu0);

  for (int i=1; i<=E.rows(); i++)
  {
    Complex kz2 = E(i)/k0/k0;
    Complex kz = sqrt(kz2);

    if (imag(kz) > 0)
      kz = -kz;

    if (abs(imag(kz)) < abs(real(kz)))
      if (real(kz) < 0)
        kz = -kz;

    cVector Ex(MN,fortranArray);
    cVector Ey(MN,fortranArray);
    cVector Hx(MN,fortranArray);
    cVector Hy(MN,fortranArray);

    Ex = eig  (r1,i);
    Ey = eig  (r2,i);
    Hx = eig_H(r1,i)/kz/k0*Y0;
    Hy = eig_H(r2,i)/kz/k0*Y0;

    BlochSectionMode* newmode 
      = new BlochSectionMode(TE_TM, kz, this, Ex, Ey, Hx, Hy);

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
// BlochSection2D::order
//
/////////////////////////////////////////////////////////////////////////////

int BlochSection2D::order(Polarisation pol, int Mx, int My) const
{
  py_error("'order' not defined for non-uniform BlochSections.");
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochSection2D::set_theta_phi
//
/////////////////////////////////////////////////////////////////////////////

void BlochSection2D::set_theta_phi(Real theta, Real phi) const
{  
  py_error("'set_theta_phi' not defined for non-uniform BlochSections.");
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformBlochSection::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void UniformBlochSection::find_modes()
{
  // Set constants.

  const Complex k = 2*pi/global.lambda*core->n();
  const Complex Y = sqrt(core->eps()/core->mu());
  
  int M = global_blochsection.Mx;
  int N = global_blochsection.My;

  int m_ = 2*M+1;  
  int n_ = 2*N+1;

  const int MN = m_*n_;

  Complex alpha0 = global_blochsection.alpha0;
  Complex  beta0 = global_blochsection.beta0;
  
  // Check values.

  if (real(global.lambda) == 0)
  {
    py_error("Error: wavelength not set.");
    return;
  }

  // Clear modeset.

  for (unsigned int i=0; i<modeset.size(); i++)
    delete modeset[i];
  modeset.clear();

  // Add modes.

  for (Real m=-M; m<=M; m+=1.0)
    for (Real n=-N; n<=N; n+=1.0)
    {
      int i1 = int((m+M+1) + (n+N)*(2*M+1));

      Complex alpha  = alpha0 + m*2.*pi/W;
      Complex  beta  =  beta0 + n*2.*pi/H;

      if ((abs(imag(alpha)) > 1e-6) || (abs(imag(beta)) > 1e-6))
        py_print("Warning:formulas not checked for complex alpha and beta.");

      Complex phi = atan2(real(beta),real(alpha));

      Complex kz = sqrt(k*k - alpha*alpha - beta*beta);

      if (imag(kz) > 0)
        kz = -kz;

      if (abs(imag(kz)) < abs(real(kz)))
        if (real(kz) < 0)
          kz = -kz;
      
      // Add TE mode.

      cVector Ex(MN,fortranArray);
      cVector Ey(MN,fortranArray);
      cVector Hx(MN,fortranArray);
      cVector Hy(MN,fortranArray);

      Ex(i1) =      -sin(phi);
      Ey(i1) =       cos(phi);
      Hx(i1) = - Y * cos(phi) * kz / k;
      Hy(i1) = - Y * sin(phi) * kz / k;

      UniformBlochSectionMode *newmode
        = new UniformBlochSectionMode(TE, kz, this, int(m), int(n), 
                                      Ex, Ey, Hx, Hy);

      newmode->normalise();
      modeset.push_back(newmode);

      // Add TM mode.

      Ex(i1) =       cos(phi) * kz / k;
      Ey(i1) =       sin(phi) * kz / k;
      Hx(i1) = - Y * sin(phi);
      Hy(i1) =   Y * cos(phi);
 
      newmode
        = new UniformBlochSectionMode(TM, kz, this, int(m), int(n), 
                                      Ex, Ey, Hx, Hy);

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
// UniformBlochSection::order
//
//  Convert polarisation and orders into the correct modal index
//
/////////////////////////////////////////////////////////////////////////////

int UniformBlochSection::order(Polarisation pol, int Mx, int My) const
{
  for (int i=0; i<modeset.size(); i++)
  {
    UniformBlochSectionMode* m 
      = dynamic_cast<UniformBlochSectionMode*>(modeset[i]);
    
    if ( (m->pol == pol) && (m->get_Mx() == Mx) && (m->get_My() == My) )
      return i;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// UniformBlochSection::set_theta_phi
//
/////////////////////////////////////////////////////////////////////////////

void UniformBlochSection::set_theta_phi(Real theta, Real phi) const
{
  Complex k = 2.*pi / global.lambda * core->n();

  global_blochsection.alpha0 = k * sin(theta) * cos(phi);
  global_blochsection.beta0  = k * sin(theta) * sin(phi);;  
}



