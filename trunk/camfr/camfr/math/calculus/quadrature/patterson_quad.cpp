
/////////////////////////////////////////////////////////////////////////////
//
// File:     patterson_quad.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000320
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "patterson.h"
#include "patterson_quad.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_sub
//
//   Helper routine for patterson_quad.
//   Integrates subinterval, but uses a relaxed convergence criterion based
//   on the estimated value of the integral over the entire interval.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad_sub(RealFunction& f, Real a, Real b, Real eps,
                        Real result_estimate, unsigned int max_k=8)
{
  bool error;
  Real abs_error;
  Real result = patterson(f, a, b, eps, &error, max_k, &abs_error);

  if ( (error == false) || (abs(abs_error) <= abs(result_estimate * eps)) )
    return result;
  else
    return patterson_quad_sub(f, a, (a+b)/2., eps, result_estimate, max_k)
         + patterson_quad_sub(f, (a+b)/2., b, eps, result_estimate, max_k);
}



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad(RealFunction& f, Real a, Real b,
                    Real eps, unsigned int max_k=8)
{
  // Try patterson on the entire interval.

  bool error;
  Real result = patterson(f, a, b, eps, &error, max_k);

  if (error == false)
    return result;

  // Do adaptive subdivision of interval.
  
  return patterson_quad_sub(f,  a, (a+b)/2., eps, result, max_k)
       + patterson_quad_sub(f, (a+b)/2., b,  eps, result, max_k);
}



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_non_adapt_sub
//
//   Helper routine for patterson_quad_non_adapt.
//   Integrates over a range of subintervals and uses a relaxed convergence
//   criterion based on the estimated value of the integral over the entire
//   interval.
//   To minimise the amount of possible work thrown away, the interval order
//   can be reversed based on the value of 'reverse' set by the observed
//   trend.
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad_non_adapt_sub
  (RealFunction& f, Real a, Real b, Real eps, Real result_estimate,
   unsigned int max_k, unsigned int N, unsigned int N_start,
   unsigned int N_stop, unsigned int* i_problem, bool* reverse)
{
  const Real delta = (b-a)/N;
  *i_problem = 0;
  Real result = 0;
  
  for (unsigned int ii=N_start; ii<=N_stop; ii++)
  {
    unsigned int i = (*reverse) ? N_start+N_stop-ii : ii;
    
    bool error;
    Real abs_error; 
    Real result_i =
      patterson(f, a+(i-1)*delta, a+i*delta, eps, &error, max_k, &abs_error);

    if ( (error==false) || (abs(abs_error) <= abs(result_estimate * eps)) )
      result += result_i;
    else
    {
      *i_problem = i;
      *reverse = (i == 2*(i/2));
      
      return 0.0;
    }
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// patterson_quad_non_adapt
//
/////////////////////////////////////////////////////////////////////////////

Real patterson_quad_non_adapt(RealFunction& f, Real a, Real b,
                              Real eps, unsigned int max_k=8)
{
  // Try patterson on the entire interval.

  bool error;
  Real result_estimate = patterson(f, a, b, eps, &error, max_k);

  if (error == false)
    return result_estimate;

  // Do non-adaptive subdivision of interval in N subintervals.

  const unsigned int max_N = 4096;
  unsigned int i_problem = 1;
  bool reverse = false;
  Real result = 0.0;

  for (unsigned int N=2; N<=max_N; N*=2)
  {
    result = 0.0;

    // Update i_problem to reflect double number of subdivisions.

    i_problem = 2*i_problem - 1;
    unsigned int i_problem_bak = i_problem;


    
    // Try problematic interval first, before starting work on other
    // intervals, work that might be discarded.

    unsigned int N_start = i_problem_bak;
    unsigned int N_stop  = i_problem_bak + 1;
    
    Real result_sub = patterson_quad_non_adapt_sub
      (f,a,b,eps, result_estimate, max_k, N,N_start,N_stop,
       &i_problem, &reverse);

    if (i_problem == 0) // No problems occured.
      result += result_sub;
    else // Discard rest of loop and double N.
      continue; 


    
    // Intervals to the left of problematic interval.

    N_start = 1;
    N_stop  = i_problem_bak - 1;
    
    result_sub = patterson_quad_non_adapt_sub
      (f,a,b,eps, result_estimate, max_k, N,N_start,N_stop,
       &i_problem, &reverse);
    
    if (i_problem == 0) // No problems occured.
      result += result_sub;
    else // Discard rest of loop and double N.
      continue;


    
    // Intervals to the right of problematic interval.

    N_start = i_problem_bak + 2;
    N_stop  = N;
    
    result_sub = patterson_quad_non_adapt_sub
      (f,a,b,eps, result_estimate, max_k, N,N_start,N_stop,
       &i_problem, &reverse);
    
    if (i_problem == 0) // No problems occured.
      result += result_sub;
    else // Discard rest of loop and double N.
      continue;


    
    // If we got here, we converged on all subintervals.
    
    return result;
  }

  cout << "Warning: maximum number of subdivision reached "
       << "in patterson_quad_non_adapt." << endl;

  return result;
}
