
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

//#include "arscomp.h"
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

SectionDisp::SectionDisp(Stack& _left, Stack& _right, Real _lambda, int _M, 
                         bool sym)
  : left(&_left), right(&_right), lambda(_lambda), M(_M), symmetric(sym) 
{
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

  // Determine minimum refractive index.

  vector<Material*> materials = left->get_materials();
  vector<Material*> right_mat = right->get_materials();
  materials.insert(materials.end(), right_mat.begin(), right_mat.end());

  min_eps_mu = materials[0]->eps_mu();
  
  for (unsigned int i=1; i<materials.size(); i++)
  {
    Complex eps_mu = materials[i]->eps_mu();
    
    if (real(eps_mu) < real(min_eps_mu))
      min_eps_mu = eps_mu;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::operator()
//
//   kt = k in the region with minimum refractive index.
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::operator()(const Complex& kt)
{
  counter++;
  
  global.lambda = lambda;
  global.orthogonal = false;
  global.polarisation = TE_TM;

  Complex old_beta = global_slab.beta;

  const Complex C = pow(2*pi/lambda, 2) / (eps0 * mu0);
  Complex beta = sqrt(C*min_eps_mu - kt*kt);

  if (real(beta) < 0)
    beta = -beta;
  
  if (abs(imag(beta)) < 1e-12)
    if (imag(beta) > 0)
      beta = -beta;

  global_slab.beta = beta;

  int old_N = global.N;
  global.N = M;

  Complex res = (global.eigen_calc==lapack) ? calc_lapack () : calc_arnoldi();

  global.N = old_N;
  global_slab.beta = old_beta;

  return res;
}



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::calc_lapack()
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::calc_lapack2()
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

  for (int i=1; i<=M; i++)
    Q(i,i) -= 1.0;
  
  cVector e(M,fortranArray);
  if (global.stability == normal)
    e.reference(eigenvalues(Q));
  else
    e.reference(eigenvalues_x(Q));
  
  // Return product of eigenvalues (determinant).

  Complex product = 1.0;

  for (int i=1; i<=M; i++)
    product *= e(i);
  
  return product;
}


Complex SectionDisp::calc_lapack()
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
  
  // Return minimum distance of eigenvalues to 1.

  int min_index = 1;

  for (int i=2; i<=M; i++)
    if (abs(e(i) - 1.0) < abs(e(min_index) - 1.0))
      min_index = i;

  return e(min_index) - 1.0;
}



/////////////////////////////////////////////////////////////////////////////
//
// Multiplier
//
//   Auxiliary class for efficient calculation of the product of a vector
//   with the cavity matrix.
//
/////////////////////////////////////////////////////////////////////////////

class Multiplier
{
  public:

    Multiplier(const cMatrix& L_, const cMatrix& R_) 
      : L(&L_), R(&R_), Q(NULL), counter(0) {}

    ~Multiplier() {delete Q;}
   
    void mult(Complex* in, Complex* out) 
    {
      counter++;
      
      cVector i(in,  global.N, blitz::neverDeleteData, fortranArray);
      cVector o(out, global.N, blitz::neverDeleteData, fortranArray); 

      if (counter == global.N)
      {
        Q = new cMatrix(global.N, global.N, fortranArray);
        Q->reference(multiply(*L, *R));
        for (int k=1; k<=global.N; k++)
          (*Q)(k,k) -= 1.0;
      }

      if (counter < global.N)
        o = multiply(*L, *R, i) - i;
      else
        o = multiply(*Q, i);
    }

  protected:

    int counter;
    
    const cMatrix* L;
    const cMatrix* R;

    cMatrix* Q;
};



/////////////////////////////////////////////////////////////////////////////
//
// SectionDisp::calc_arnoldi()
//
/////////////////////////////////////////////////////////////////////////////

Complex SectionDisp::calc_arnoldi()
{
/*
  left->calcRT();
  if (! symmetric)
    right->calcRT();

  Multiplier m (              left->as_multi()->get_R12(),
                 symmetric ?  left->as_multi()->get_R12()
                           : right->as_multi()->get_R12());

  int nev = 1;
  int ncv = 6;
  double tol = 1e-3;
  int max_iter = 1000;

  ARCompStdEig<Real, Multiplier> prob
    (M, nev, &m, &Multiplier::mult, "SM", ncv, tol, max_iter);

  prob.FindEigenvalues();

  if (prob.ConvergedEigenvalues() == 0)
  {
    cout << "Warning: Arnoldi solver did not converge for beta "
         << beta << endl << "Using LAPACK instead. " << endl;
    
    return calc_lapack(beta);
  }

  return prob.Eigenvalue(0);
*/
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

  params.push_back(min_eps_mu);
  
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

  min_eps_mu = params[params_index];

   left->freeRT();
  right->freeRT();
}
