
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectiondisp.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020125
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "section.h"
#include "sectiondisp.h"
#include "../slab/generalslab.h"
#include "../slab/slabmatrixcache.h"
#include "../../math/linalg/linalg.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::SectionDisp()
//
/////////////////////////////////////////////////////////////////////////////

SectionDisp::SectionDisp(Stack& _left, Stack& _right, const Complex& _lambda, 
                         int _M, bool sym)
  : left(&_left), right(&_right), lambda(_lambda), M(_M), symmetric(sym) 
{
  // TODO: Rework.

  // Determine left slabs.

  const vector<Chunk>* l 
    = dynamic_cast<StackImpl*>(left->get_flat_sc())->get_chunks();

  for (unsigned int i=0; i<l->size(); i++)
  {
    left_slabs.push_back(dynamic_cast<Slab*>((*l)[i].sc->get_inc()));
    left_slabs.push_back(dynamic_cast<Slab*>((*l)[i].sc->get_ext()));
  }

  remove_copies(&left_slabs);

  // Determine right slabs.

  const vector<Chunk>* r 
    = dynamic_cast<StackImpl*>(right->get_flat_sc())->get_chunks();

  for (unsigned int i=0; i<r->size(); i++)
  {
    right_slabs.push_back(dynamic_cast<Slab*>((*r)[i].sc->get_inc()));
    right_slabs.push_back(dynamic_cast<Slab*>((*r)[i].sc->get_ext()));
  } 

  remove_copies(&right_slabs);

  // Get pointers to top-level (not flat_sc) thicknesses.

  left_d  = dynamic_cast<StackImpl*>( left->get_sc())->get_thicknesses();
  right_d = dynamic_cast<StackImpl*>(right->get_sc())->get_thicknesses();

  // Determine slabs.

  for (int i=l->size()-1; i>=0; i--)
  {
    Complex d_i = (*l)[i].d;

    if (abs(d_i) > 1e-6)
    {
      d.push_back(d_i);
      slabs.push_back(dynamic_cast<Slab*>((*l)[i].sc->get_ext()));
    }
  }

  for (int i=0; i<r->size(); i++)
  {
    if (i==0 && (*r)[i].sc->get_ext()==slabs.back())
    {
      d.back() += (*r)[i].d;
      continue;
    }
    
    Complex d_i = (*r)[i].d;

    if (abs(d_i) > 1e-6)
    {
      d.push_back(d_i);    
      slabs.push_back(dynamic_cast<Slab*>((*r)[i].sc->get_ext()));
    }
  }

  // Determine refractive index of material which is best used for kt:
  // minimum for dielectics, maximum for metals.

  vector<Material*> materials = left->get_materials();
  vector<Material*> right_mat = right->get_materials();
  materials.insert(materials.end(), right_mat.begin(), right_mat.end());

  vector<Complex> eps_mu;
  for (unsigned int i=0; i<materials.size(); i++)
    eps_mu.push_back(materials[i]->eps_mu());

  remove_copies(&eps_mu, 1e-23);
  std::sort(eps_mu.begin(), eps_mu.end(), RealSorter());
  
  Complex min_eps_mu, max_eps_mu;
  
  if (eps_mu.size() >= 2)
  {
     min_eps_mu = eps_mu[0];
     max_eps_mu = eps_mu[eps_mu.size()-1];
  }
  else
    min_eps_mu = max_eps_mu = eps_mu[0];

  //bool metallic = (real(min_eps_mu) < 0);
  bool metallic = global_slab.low_index_core;

  kt_eps_mu = metallic ? max_eps_mu : min_eps_mu;
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::operator()
//
//   kt = k in the region with minimum refractive index, or in the region 
//   with maximum refractive index if metallic structures are involved
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::operator()(const Complex& kt2)
{
    counter++;

    global.lambda = lambda;
    global.polarisation = TE_TM;

    bool old_orthogonal = global.orthogonal;
    global.orthogonal = false;

    Complex old_beta = global.slab_ky;

    const Complex C = pow(2*pi/lambda, 2) / (eps0 * mu0);
    Complex beta = sqrt(C*kt_eps_mu - kt2);

    if (real(beta) < 0)
        beta = -beta;
  
    if (abs(imag(beta)) < 1e-12)
        if (imag(beta) > 0)
            beta = -beta;

    global.slab_ky = beta;

    int old_N = global.N;
    global.N = M;
   
    Complex res = (global.eigen_calc==lapack) ? calc_split () : calc_global();
   
    global.N = old_N;
    global.slab_ky = old_beta;
    global.orthogonal = old_orthogonal;
   
    return res;
}



/////////////////////////////////////////////////////////////////////////////
//
// Index mapping into dense (ku=-1) or band (ku=upper diagonals) storage.
//
/////////////////////////////////////////////////////////////////////////////

struct IndexMap
{
  int ku;
  IndexMap(int ku_) : ku(ku_) {}
  int operator()(int i, int j) {return (ku == -1) ? i : ku+1+i-j;}  
};



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::calc_global()
//
//   Handles all layers of the stack at the same time.
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::calc_global()
{
  // Allocate matrices.

  const int K = slabs.size();
  int kl, ku; kl = ku = 3*M-1; // lower and upper diagonals.
  bool band_storage = (K > 3 - .5/M); // only quicker if 6M-1 < 2MK

  // Since the new version uses eigenvalues instead of determinants, we
  // disable band storage, as we don't have a eigensolver for band matrices.

  band_storage = false;
  
  IndexMap f(band_storage ? ku : -1);

  cMatrix Q(band_storage ? kl+1+ku : 2*M*K, 2*M*K, fortranArray);
  Q = 0.0;

  cVector prop_I(M,fortranArray), prop_II(M,fortranArray);

  // Reflection from left wall.

  Complex R0 = (global_section.leftwall == E_wall) ? -1.0 : 1.0;

  slabs[0]->find_modes();

  for (int i=1; i<=M; i++)
  {
    Q(f(i,  i), i  ) = 1.0;
    Q(f(i,M+i), M+i) = -R0 * exp(-I*slabs[0]->get_mode(i)->get_kz() * d[0]);
  }

  // Interfaces.

  for (unsigned int k=1; k<K; k++)
  {
    // Calculate propagation vectors.

    slabs[k]->find_modes();

    for (int i=1; i<=M; i++)
    {
      prop_I (i) = exp(-I*slabs[k-1]->get_mode(i)->get_kz() * d[k-1]);
      prop_II(i) = exp(-I*slabs[k  ]->get_mode(i)->get_kz() * d[k  ]);      
    }
   
    // Calculate overlap matrices.

    cMatrix O_I_I (M,M,fortranArray); cMatrix O_II_II(M,M,fortranArray);
    cMatrix O_I_II(M,M,fortranArray); cMatrix O_II_I (M,M,fortranArray);

    slabs[k-1]->calc_overlap_matrices(slabs[k],
                                      &O_I_II, &O_II_I, &O_I_I, &O_II_II);

    // Index of block's top left row and column (in dense storage).

    int r = (k-1)*2*M + M;
    int c = (k-1)*2*M;

    // Fill matrices.

    for (int i=1; i<=M; i++)
      for (int j=1; j<=M; j++)
      {
        Q(f(  r+i,    c+j),     c+j) =  O_I_II(j,i)  * prop_I(j);
        Q(f(M+r+i,    c+j),     c+j) =  O_II_I(i,j)  * prop_I(j);

        Q(f(  r+i,  M+c+j),   M+c+j) =  O_I_II(j,i);
        Q(f(M+r+i,  M+c+j),   M+c+j) = -O_II_I(i,j);

        Q(f(  r+i,2*M+c+j), 2*M+c+j) = -O_II_II(j,i);
        Q(f(M+r+i,2*M+c+j), 2*M+c+j) = -O_II_II(i,j);

        Q(f(  r+i,3*M+c+j), 3*M+c+j) = -O_II_II(j,i) * prop_II(j);
        Q(f(M+r+i,3*M+c+j), 3*M+c+j) =  O_II_II(i,j) * prop_II(j);
      } 
  }

  // Reflection from right wall.

  Complex Rn = (global_section.rightwall == E_wall) ? -1.0 : 1.0;

  int r = 2*M*K-  M;
  int c = 2*M*K-2*M;

  for (int i=1; i<=M; i++)
  {
    Q(f(r+i,  c+i),  c+i)=-Rn*exp(-I*slabs[K-1]->get_mode(i)->get_kz()*d[K-1]);
    Q(f(r+i,M+c+i),M+c+i)= 1.0;
  }

  cVector e(Q.rows(),fortranArray);
  if (global.stability == normal)
    e.reference(eigenvalues(Q));
  else
    e.reference(eigenvalues_x(Q));
  
  // Return product of K best eigenvalues (i.e. closest to 0), 
  // but don't count the lowest two. (These are parasitic.)

  int K_ = 3;

  Complex product = 1.0;
  vector<unsigned int> min_indices;

  for (unsigned int k=0; k<K_; k++)
  {
    int min_index = 1;
    Real min_distance = 1e200;
    for (int i=1; i<=Q.rows(); i++)
    {
      if ( (abs(e(i)) < min_distance) &&
           (std::find(min_indices.begin(),min_indices.end(),i) 
            == min_indices.end()))
      {
        min_index = i;
        min_distance = abs(e(i));
      } 
    }

    // Don't count eigenvalues twice.

    if (std::find(min_indices.begin(),min_indices.end(),min_index) 
         == min_indices.end())
    {
      product *= e(min_index);
      min_indices.push_back(min_index);
    }
  }

  // Only return highest of the three.

  int final_index = 0;
  Real max_norm = 0;
 
  for (int i=0; i< min_indices.size(); i++)
  {
    if (abs(e(min_indices[i])) > max_norm)
    {
      final_index = min_indices[i];
      max_norm = abs(e(min_indices[i]));
    }
  }

  return e(final_index);

  // Old version: return the determinant.

  // This is prone to overflow, especially for a large number of layers.
  // TODO: try Arnoldi eigenvalue solver for this, or Rayleigh scheme from
  // TOMS.

  if (band_storage)
    return determinant_band(Q, 2*K*M, kl, ku);
  else
    return determinant(Q);
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::calc_split()
//
//   Splits the stack in a left and right stack.
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::calc_split()
{
  // Calculate eigenvectors.
 
  left->calcRT();
  if (! symmetric)
    right->calcRT();
  
  cMatrix Q(M,M,fortranArray);
  if (! symmetric)
    Q.reference(multiply( left->as_multi()->get_R12(), 
                         right->as_multi()->get_R12()));
  else
    Q.reference(multiply( left->as_multi()->get_R12(), 
                          left->as_multi()->get_R12()));
  
  cVector e(M,fortranArray);
  if (global.stability == normal)
    e.reference(eigenvalues(Q));
  else
    e.reference(eigenvalues_x(Q));
  
  // Return product of K best eigenvalues (i.e. closest to 1).
  // Don't take zero eigenvalues into account as these are not numerically
  // relevant.

  int K = 5;

  Complex product = 1.0;
  vector<unsigned int> min_indices;

  for (unsigned int k=0; k<K; k++)
  {
    int min_index = 1;
    Real min_distance = 1e200;
    for (int i=1; i<=M; i++)
    {
      if ( (abs(e(i) - 1.0) < min_distance) &&
           (std::find(min_indices.begin(),min_indices.end(),i) 
            == min_indices.end()) &&
            (abs(e(i)) > 1e-5) )
      {
        min_index = i;
        min_distance = abs(e(i) - 1.0);
      } 
    }

    // Don't count eigenvalues twice.

    if (std::find(min_indices.begin(),min_indices.end(),min_index) 
         == min_indices.end())
    {
      product *= e(min_index) - 1.0;
      min_indices.push_back(min_index);
    }
  }

  return product;
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::calc_T()
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::calc_T()
{
  // Calculate eigenvectors.
 
  left->calcRT();
  if (! symmetric)
    right->calcRT();

  cMatrix Q(M,M,fortranArray);
  Q = left->as_multi()->get_T12();
  
  cVector e(M,fortranArray);
  if (global.stability == normal)
    e.reference(eigenvalues(Q));
  else
    e.reference(eigenvalues_x(Q));  
  
  // Return product of K best eigenvalues (i.e. largest).
  // Don't take zero eigenvalues into account as these are not numerically
  // relevant.

  int K = 1; //5;

  Complex product = 1.0;
  vector<unsigned int> max_indices;

  for (unsigned int k=0; k<K; k++)
  {
    int max_index = 1;
    Real max_distance = 1e-200;
    for (int i=1; i<=M; i++)
    {
      if ( (abs(e(i)) > max_distance) &&
           (std::find(max_indices.begin(),max_indices.end(),i) 
            == max_indices.end()) )
      {
        max_index = i;
        max_distance = abs(e(i));
      } 
    }

    // Don't count eigenvalues twice.

    if (std::find(max_indices.begin(),max_indices.end(),max_index) 
         == max_indices.end())
    {
      product *= 1.0/e(max_index);
      max_indices.push_back(max_index);
    }
  }

  return product;
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::get_params()
//
///////////////////////////////////////////////////////////////////////////// 

vector<Complex> SectionDisp::get_params() const
{
  vector<Complex> params;

  for (unsigned int i=0; i<left_slabs.size(); i++)
  {
    vector<Complex> params_i = left_slabs[i]->get_params();
    params.push_back(params_i.size());
    params.insert(params.end(), params_i.begin(), params_i.end());
  }

  for (unsigned int i=0; i<right_slabs.size(); i++)
  {
    vector<Complex> params_i = right_slabs[i]->get_params();
    params.push_back(params_i.size());
    params.insert(params.end(), params_i.begin(), params_i.end());
  }

  for (unsigned int i=0; i<left_d.size(); i++)
    params.push_back(*left_d[i]);

  for (unsigned int i=0; i<right_d.size(); i++)
    params.push_back(*right_d[i]);

  params.push_back(lambda);

  params.push_back(kt_eps_mu);
  
  return params;
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::set_params()
//
/////////////////////////////////////////////////////////////////////////////

void SectionDisp::set_params(const vector<Complex>& params)
{
  unsigned int params_index = 0;

  for (unsigned int i=0; i<left_slabs.size(); i++)
  {
    unsigned int j_max = (unsigned int)(real(params[params_index++]));

    vector<Complex> params_i;
    for (unsigned int j=0; j<j_max; j++)
      params_i.push_back(params[params_index++]);
    
    left_slabs[i]->set_params(params_i);
  }

  for (unsigned int i=0; i<right_slabs.size(); i++)
  {
    unsigned int j_max = (unsigned int)(real(params[params_index++]));

    vector<Complex> params_i;
    for (unsigned int j=0; j<j_max; j++)
      params_i.push_back(params[params_index++]);
    
    right_slabs[i]->set_params(params_i);
  }

  for (unsigned int i=0; i<left_d.size(); i++)
    *left_d[i]  = params[params_index++];

  for (unsigned int i=0; i<right_d.size(); i++)
    *right_d[i] = params[params_index++];
  
  lambda = real(params[params_index++]);

  kt_eps_mu = params[params_index];

   left->freeRT();
  right->freeRT();
}
