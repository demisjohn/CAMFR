
/////////////////////////////////////////////////////////////////////////////
//
// File:     root.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990513
// Version:  1.1
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include "root.h"

/*
 ************************************************************************
 *
 *			  Numerical Math Package : http://pobox.com/~oleg/ftp
 *
 *			    Brent's root finder
 *	       obtains a zero of a function of one variable
 *        with accuracy 4*machine_eps()*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 * The function makes use of a bissection procedure combined with
 * a linear or quadratic inverse interpolation.
 * At each step the code operates three abscissae - a, b, and c:
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even an earlier approximation such that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c encompass
 *		   the root
 * Given these abscissae, the code computes two new approximations, one by the 
 * bissection procedure and the other one from interpolation (if a,b, and c
 * are all different the quadratic interpolation is used, linear otherwise).
 * If the approximation obtained by the interpolation looks
 * reasonable (i.e. falls within the current interval [b,c], not too close
 * to the end points of the interval), the point is accepted as a new
 * approximation to the root. Otherwise, the result of the bissection is used.
 * Therefore, the range of uncertainty is guaranteed to tighten at 
 * least by a factor of 1.6
 *
 * $Id: root.cpp,v 1.2 2002-02-22 20:45:58 pbienst Exp $
 *
 ************************************************************************
 */

Real brent_root(Function1D<Real>& f, Real ax, Real bx, Real eps)
{
  const int MAXITER = 100;
  static const Real mach_eps = machine_eps();

  if (eps < 0)
  {
    cerr << "Tolerance must be positive." << endl;
    exit (-1);
  }

  if (ax > bx)
  {
    cerr << "Left end point of the interval should be strictly less "
         << "than the right one." << endl;
    exit (-1);
  }
  
  Real  b = bx;   // last and best approximation to the root
  Real fb = f(b);
  Real  a = ax;   // last but one approximation
  Real fa = f(a);
  Real  c = a;    // last but one or even an earlier approx
  Real fc = fa;   // (see the condition above)

  if (fa*fb > 0)
  {
    cerr << "Error: no certain root in interval ["
         << a << "," << b <<"]" << endl;
    exit (-1);
  }
  
  for (int i=0; i<MAXITER; i++)
  {
    
    const Real prev_step = b-a;

    // swap data to make b the best approximation found so far
    
    if ( abs(fc) < abs(fb) )
    {                            
      a = b;  b = c;  c = a;
      fa=fb;  fb=fc;  fc=fa;
    }

    const Real eps_act = 2*mach_eps*abs(b) + eps/2;
    Real new_step = (c-b)/2;

    if ( (abs(new_step) <= eps_act) || (abs(fb) <= eps_act) )
      return b;

    // figure out if interpolation can be tried
    
    if ( (abs(prev_step) >= eps_act) && (abs(fa) > abs(fb)) )  
    {                           
      Real p,q;             
      const Real cb = c-b;

      if ( a==c ) // do linear interpolation
      {					        
        register const Real t1 = fb/fa;
        p = cb*t1;
        q = 1.0 - t1;
      }
      else			// do quadratic inverse interpolation
      {
        register const Real t1=fb/fc, t2=fb/fa;
        q = fa/fc;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }

      // Formulas above computed new_step = p/q with the wrong sign
      // (on purpose). Correct this, but in such a way so that p
      // would be positive.
      
      if ( p > 0 )
        q = -q;
      else
        p = -p;

      // If b+p/q falls in [b,c] and isn't too large it is accepted.
      // If p/q is too large then the bissection procedure can
      // reduce [b,c] to a larger extent.
      
      if ( (2*p < (1.5*cb*q-abs(eps_act*q))) && (2*p < abs(prev_step*q)) )
        new_step = p/q;
    }

    // Adjust the step to be not less than the tolerance
    
    if ( abs(new_step) < eps_act )	     
      new_step =  (new_step > 0) ? eps_act : -eps_act;

    // Update values for next iteration
    
    a = b;  fa = fb;
    b += new_step;  fb = f(b);

    // Adjust c for it to have the sign opposite to that of b
    
    if ( ((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0)) )
    {                 			           
      c = a;  fc = fa;
    }
    
  }

  cout << "Warning: Brent solver did not reach requested accuracy " << eps 
       << " in interval [" << ax << "," << bx << "]. " << endl;

  return b;
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_root
//
//  same as above, but can handle different intervals at the same time
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_root(Function1D<Real>& f, vector<Real>& Ax,
                        vector<Real>& Bx, Real eps)
{
  if (Ax.size() != Bx.size())
  {
    cerr << "Error: Ax and Bx must be the same size." << endl;
    exit (-1);
  }

  vector<Real> zeros;
  
  for (unsigned int i=0; i<Ax.size(); i++)
    zeros.push_back(brent_root(f, Ax[i], Bx[i], eps));
  
  return zeros;
}

    
    
/////////////////////////////////////////////////////////////////////////////
//
// bracket_all_roots
//
/////////////////////////////////////////////////////////////////////////////

void bracket_all_roots(Function1D<Real>& f, Real ax, Real bx,
                       vector<Real>& Ax, vector<Real>& Bx, Real dx,
                       int sec_level)
{ 
  Ax.clear(); Bx.clear();
    
  Real fx, previous_fx;
  Real fx_coarse, previous_fx_coarse;

  fx = previous_fx = fx_coarse = previous_fx_coarse = f(ax);

  long int iters = 0;
  int coarse_zeros = 0;
  
  for (Real x=ax+dx; x<=bx; x+=dx, iters++)
  {
    fx = f(x);
    
    if (fx * previous_fx <= 0)
    {
      Ax.push_back(x-dx);
      Bx.push_back(x);
    }
    previous_fx = fx;

    if (sec_level > 0)
    {
      if ( ((iters % int(pow(2.0,sec_level)))==0) || (x+dx>bx) )
      {
        fx_coarse = fx;
        if (fx_coarse * previous_fx_coarse <= 0)
          coarse_zeros++;
        previous_fx_coarse = fx_coarse;
      }
    }
    
  }

  if ( (sec_level > 0) && (int(Ax.size()) != coarse_zeros) )
  {
    cout << "Warning: number of zeros with coarse grid different than "
         << "with fine grid. " << endl;
    cout << "Step   dx: " << Ax.size() << endl;
    cout << "Step " << pow(2.0,sec_level) << "*dx: " << coarse_zeros << endl;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// bracket_N_roots
//
/////////////////////////////////////////////////////////////////////////////

void bracket_N_roots(Function1D<Real>& f, Real ax, int N,
                     vector<Real>& Ax, vector<Real>& Bx,
                     Real dx, int sec_level)
{
  Ax.clear(); Bx.clear();

  Real fx, previous_fx;
  Real fx_coarse, previous_fx_coarse;

  fx = previous_fx = fx_coarse = previous_fx_coarse = f(ax);

  const long int MAXITER = 10000000L;
  long int iters         = 0;
       int fine_zeros    = 0;
       int coarse_zeros  = 0;       
  
  for (Real x=ax+dx; fine_zeros<N; x+=dx, iters++)
  { 
    fx = f(x);
    
    if (fx * previous_fx <= 0)
    {      
      Ax.push_back(x-dx);
      Bx.push_back(x);
      fine_zeros++;
    }
    
    previous_fx = fx;

    if (sec_level > 0)
    {
      if (fine_zeros==N)
        coarse_zeros++; // Avoid problems at the end of the interval.
      else if ( (iters % int(pow(2.0,sec_level)))==0 )
      {
        fx_coarse = fx;
        if (fx_coarse * previous_fx_coarse <= 0)
          coarse_zeros++;
        previous_fx_coarse = fx_coarse;
      }
    }
    
  }

  if (iters >= MAXITER-1)
    cout << "Warning: maximum iterations " << MAXITER << " reached." << endl;

  if ( (sec_level > 0) && (int(Ax.size()) != coarse_zeros) )
  {
    cout << "Warning: number of zeros with coarse grid different than "
         << "with fine grid. " << endl;
    cout << "Step   dx: " << Ax.size() << endl;
    cout << "Step " << pow(2.0,sec_level) << "*dx: " << coarse_zeros << endl;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_all_roots
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_all_roots(Function1D<Real>& f, Real ax, Real bx, Real dx,
                             Real eps, int sec_level)
{
  vector<Real> Ax, Bx; 
  bracket_all_roots(f, ax, bx, Ax, Bx, dx, sec_level);
  return brent_root(f, Ax, Bx, eps);
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_N_roots
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_N_roots(Function1D<Real>& f, Real ax, int N, Real dx,
                           Real eps=1e-13, int sec_level=0)
{
  vector<Real> Ax, Bx;
  bracket_N_roots(f, ax, N, Ax, Bx, dx, sec_level);
  return brent_root(f, Ax, Bx, eps); 
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_refine_roots
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_refine_roots(Function1D<Real>& f, vector<Real>& x,
                                Real delta_x, Real dx,
                                Real eps, int sec_level)
{
  sort(x.begin(), x.end());
  
  vector<Real> results;
  
  for (unsigned int i=0; i<x.size(); i++)
  {
    Real ax = (1-delta_x)*x[i];
    
    while ( (i<x.size()-1) && ( (1+delta_x)*x[i] >= (1-delta_x)*x[i+1] ))
       i++;
    
    Real bx = (1+delta_x)*x[i];
    
    vector<Real> x_i
      = brent_all_roots(f, ax, bx, dx, eps, sec_level);

    for (unsigned int j=0; j<x_i.size(); j++)
      results.push_back(x_i[j]);
  }
  
  return results;
}

