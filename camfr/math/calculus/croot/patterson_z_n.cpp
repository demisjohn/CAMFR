
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson_z_n.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "../function.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// patterson_z_n
//
//   This quadrature program uses formulae due to T. N. L. Patterson,
//   Mathematics of computation, Volume 22, 1968, pages 847-856, as
//   modified by F. T. Krogh and W. V. Snyder, ACM Transactions on
//   Mathematical Software 17, 4 (December 1991) pp 457-461.  It is a
//   functional replacement for Algorithm 468, T. N. L. Patterson,
//   Communications of the ACM 16, 11 (November 1973) 694-699.
//
//   The complex curve is taken to be the line segment between a and b
//   and is parameterised by an independent variable t running from 0 to 1.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> patterson_z_n(ComplexFunction& f,
                              const Complex& a, const Complex& b, int M,
                              Real eps, Real mu, bool* error_ptr,
                              unsigned int max_k,
                              vector<Complex>* abs_error)
{
  // Check if a and b are different.

  if (1. + abs(a-b) <= 1.)
  {
    vector<Complex> result;

    // TODO: - missing abs_error definition
    for (unsigned int i=0; i<=M; i++)
      result.push_back(0.0);
    
    return result;
  }
  
  // Check and coerce k.
  
  if (max_k < 2)
  {
    py_print("Warning: increasing max_k to 2.");
    max_k = 3;
  }
    
  if (max_k > 8)
  {
    py_print("Warning: restricting max_k to 8.");
    max_k = 8;
  }
  
  // Include the array 'p' with the coefficients used in these formulas. 

  #include "patterson_coeff.cpp"

  // Define constants and workspace containing previous function evaluations.

  Complex* work[18];
  for (unsigned int i=0; i<18; i++)
    work[i] = new Complex[M+1];

  Complex* fz     = new Complex[M+1];
  Complex* acum   = new Complex[M+1];
  Complex* result = new Complex[M+1];
  
  // Apply 1-point Gauss formula (midpoint rule).

  const Complex delta = b-a;
  fz[0] = 1.0 / f(a + 0.5*delta);

  //if (abs(fz[0]) > 1e3)
  //  cout << "Possible zero for " << a+0.5*delta << endl;
  
  result[0] = acum[0] = work[1][0] = fz[0];

  for (int n=1; n<=M; n++)
  {
    fz[n] = fz[n-1]*(a + 0.5*delta);
    
    result[n] = acum[n] = work[1][n] = fz[n];
    
    if (machine_eps()*abs(fz[n]) > mu)
      M = n-1;
  }
  
  // Go on to the next formulas.

  // fl, fh contain subscripts to indicate where to store integrand
  //        samples in 'work'.
  // kl, kh are bounds for indexes at which integrand samples are
  //        retrieved from 'work' in order to begin applying a quadrature
  //        formula.
  // kx     is a list of bounds of subscripts into kh and kl.  kh is
  //        indexed by k.  The set of indices from which to retrieve
  //        integrand samples is given by the set of bounds in kl and kh
  //        indexed by kx[k-1]+1 to kx[k] inclusively.

  static const int fl[] = {0, 0, 2, 3, 5, 9,12,14, 1};
  static const int fh[] = {0, 0, 2, 4, 8,16,17,17, 0};
  static const int kl[] = {0,       1, 1, 1, 1, 1, 3, 5, 9, 5, 9,12};
  static const int kh[] = {0,       1, 2, 4, 8,16, 3, 6,17, 5, 9,17};
  static const int kx[] = {0, 0,    1, 2, 3, 4, 5,       8,      11};

  int ip=1; // Index in coefficient array 'p'.
  int jh=0;

  Complex* prev_result = new Complex[M+1];
  Complex* prev_acum   = new Complex[M+1];
  
  for (int k=2; k<=max_k; k++)
  {
    for (unsigned int n=0; n<=M; n++)
    {
      prev_result[n] = result[n];
      prev_acum[n]   = acum[n];
      acum[n] = 0.0;
    }
    
    // Compute contribution to current estimate due to function
    // values used in previous formulas.
    
    for (int kk=kx[k-1]+1; kk<=kx[k]; kk++)
      for (int j=kl[kk]; j<=kh[kk]; j++)
      {
        for (unsigned int n=0; n<=M; n++)
          acum[n] += p[ip] * work[j][n];

        ip++;
      }
    
    // Compute contribution from new function values.
    
    int jl=jh+1;  
        jh=jl+jl-1;
    int j1=fl[k];
    int j2=fh[k];

    for (int j=jl; j<=jh; j++)
    {
      Complex t = p[ip++]*0.5;
      Complex f1 = 1.0 / f(a +      t *delta);
      Complex f2 = 1.0 / f(a + (1.0-t)*delta);

      //if (abs(f1) > 1e3)
      //  cout << "Possible zero for " << a+t*delta << endl;
      
      //if (abs(f2) > 1e3)
      //  cout << "Possible zero for " << a+(1.0-t)*delta << endl;
      
      fz[0] = f1+f2;
      acum[0] += p[ip]*fz[0];
      
      if (j1 <= j2)
        work[j1][0] = fz[0];
      
      for (unsigned int n=1; n<=M; n++)
      {
        f1 *= a +      t *delta;
        f2 *= a + (1.0-t)*delta;
        
        fz[n] = f1+f2;
        acum[n] += p[ip]*fz[n];

        if (j1 <= j2)
          work[j1][n] = fz[n];

        if ( (machine_eps()*abs(f1) > mu) || (machine_eps()*abs(f2) > mu) )
          M = n-1;
      }

      ip++;
      j1++;
    }

    // Accumulate result.

    bool all_converged = true;
    for (unsigned int n=0; n<=M; n++)
    {
      result[n] = acum[n] = 0.5*(acum[n] + prev_acum[n]);
      
      if (abs(result[n]-prev_result[n]) > abs(eps*result[n]))
        all_converged = false;
    }

    if (all_converged)
    {
      if (error_ptr)
        *error_ptr = false;

      goto final;
    }
  }

  // Failed to converge.

  if (error_ptr)
    *error_ptr = true;

  // Finalise.

  final:

  if (abs_error)
  {
    abs_error->clear();
    for (unsigned int n=0; n<=M; n++)
      abs_error->push_back(delta*(result[n]-prev_result[n]));
  }

  vector<Complex> final;
  for (unsigned int n=0; n<=M; n++)
    final.push_back(delta*result[n]);

  for (unsigned int i=0; i<18; i++)
    delete [] work[i];

  delete [] fz;
  delete [] acum;
  delete [] result;
  delete [] prev_result;
  delete [] prev_acum;    

  return final;
}



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_z_n_sub
//
//   Helper routine for patterson_quad_z_n.
//   Integrates subinterval, but uses a relaxed convergence criterion based
//   on the estimated value of the integral over the entire interval.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> patterson_quad_z_n_sub
  (ComplexFunction& f, const Complex& a, const Complex& b, int M,
   Real eps, Real mu, const vector<Complex>& result_estimate,
   unsigned int max_k)
{ 
  bool error;
  vector<Complex> abs_error;
  vector<Complex> result
    = patterson_z_n(f,a,b,M,eps,mu,&error,max_k,&abs_error);

  if (error == false)
    return result;

  // Check if relaxed convergence is satisfied.

  bool converged = true;

  for (unsigned int i=0; i<abs_error.size(); i++)
    // CHECK - we need to have result_estimate size >= abs_error size
    if ( abs(abs_error[i]) > abs(result_estimate[i] * eps) )
      converged = false;

  if (converged)
    return result;

  // Do adaptive subdivision of interval.

  // CHECK: M result size can be < M - so we need to call following functions with updated M otherwise it will
  // crash when estimating abs_error
  vector<Complex> result1
    = patterson_quad_z_n_sub(f, a, (a+b)/2., result.size()-1,eps,mu,result_estimate,max_k);

  vector<Complex> result2
    = patterson_quad_z_n_sub(f, (a+b)/2., b, result.size()-1,eps,mu,result_estimate,max_k);

  unsigned int new_M
    = (result1.size()>result2.size()) ? result2.size() : result1.size();

  vector<Complex> new_result;
  for (unsigned int i=0; i<new_M; i++)
    new_result.push_back(result1[i] + result2[i]);

  return new_result;
}



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_z_n
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> patterson_quad_z_n(ComplexFunction& f,
                                   const Complex& a, const Complex& b, int M,
                                   Real eps, Real mu, unsigned int max_k)
{
  // Try patterson on the entire interval.

  bool error;
  vector<Complex> result 
    = patterson_z_n(f, a, b, M, eps, mu, &error, max_k, NULL);

  if (error == false)
    return result;

  // Do adaptive subdivision of interval

  // CHECK: M result size can be < M - so we need to call following functions with updated M otherwise it will
  // crash when estimating abs_error
  vector<Complex> result1
    = patterson_quad_z_n_sub(f, a, (a+b)/2., result.size()-1, eps, mu, result, max_k);

  vector<Complex> result2
    = patterson_quad_z_n_sub(f, (a+b)/2., b, result.size()-1, eps, mu, result, max_k);

  unsigned int new_M
    = (result1.size()>result2.size()) ? result2.size() : result1.size();

  vector<Complex> new_result;
  for (unsigned int i=0; i<new_M; i++)
    new_result.push_back(result1[i] + result2[i]);

  return new_result;
}
