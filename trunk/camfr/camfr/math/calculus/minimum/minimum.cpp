
/////////////////////////////////////////////////////////////////////////////
//
// File:     minimum.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000125
// Version:  1.1
//
// Copyright (C) 1999-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include "minimum.h"
       
/*
 ************************************************************************
 *
 *			  Numerical Math Package : http://pobox.com/~oleg/ftp
 *
 *		    Brent's one-dimensional minimizer 
 *
 *	Finds a local minimum of a single argument function
 * over a given range with accuracy 3*sqrt(machine_eps())*abs(x) + eps
 *
 *	The procedure can only determine a local minimum, which coincides with
 *	the global one if and only if the function under investigation is
 *	unimodular.
 *	If a function being examined possesses no local minimum within
 *	the given interval, Fminbr returns either the left or the right end
 *	point of the interval, wherever the function value is smaller.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.202 of the Russian edition
 *
 * The function makes use of a "golden section" procedure combined with
 * a parabolic interpolation.
 * At each step the code operates three abscissae - x,v, and w.
 * 	x - the last and the best approximation to the minimum location,
 *		i.e. f(x) <= f(a) or/and f(x) <= f(b)
 * 	    (if the function f has a local minimum in (a,b), then both
 *  	     conditions are met after one or two steps).
 *	v,w are previous approximations to the location of the minimum.
 *	They may coincide with a, b, or x (although the algorithm tries
 *	to make all u, v, and w distinct). 
 * Points x, v, and w are used to construct an interpolating parabola,
 * whose minimum is regarded as a new approximation to the minimum
 * of the function, provided the parabola's minimum falls within [a,b]
 * and reduces the current interval [a,b] to a larger extent than the
 * gold section procedure does.
 * When f(x) has a positive second derivative at the point of minimum
 * (which does not coincide with a or b) the procedure converges
 * superlinearly at a rate of about 1.324
 *
 * $Id: minimum.cpp,v 1.2 2002-02-22 20:45:58 pbienst Exp $
 *
 ************************************************************************
 */

Real brent_minimum(Function1D<Real>& f, Real ax, Real bx, Real eps)
{  
  const int MAXITER = 100;
  static const Real r = (3-sqrt(5.0))/2;	// The golden section ratio.
  static const Real sqrt_mach_eps = sqrt(machine_eps());
  
  if (eps < 0)
  {
    cerr << "Tolerance must be positive." << endl;
    exit (-1);
  }

  if (ax > bx)
  {
    cout << "Left end point of the interval should be strictly less "
         << "than the right one." << endl;
    exit (-1);
  }
  
  Real a = ax, b = bx;  // Current interval.
  Real v = a + r*(b-a); // First step - always golden section.
  Real fv = f(v);
  Real x = v;				// The last and the best approximation.
  Real fx = fv;
  Real w = v;				// A previous approx to the minimum.
  Real fw = fv;

  for (int i=0; i<MAXITER; i++)
  {
    const Real range    = b-a;
    const Real midpoint = (a+b)/2;
    const Real eps_act  = sqrt_mach_eps*abs(x) + eps/3;
    
    if ( 2*abs(x-midpoint) + range <= 4*eps_act )
      return x;
    
    // Compute a new step with the golden section.
    
    Real new_step = r * ( x < midpoint ? b-x : a-x );

    // Decide on the interpolation. If x and w are distinct
    // interpolatiom may be tried.
    
    if ( abs(x-w) >= eps_act ) 
    {				
      register Real p; 	
      register Real q;
      register Real t;

      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);

      // Formulas above computed new_step = p/q with the wrong sign
      // (on purpose). Correct this, but in such a way so that p
      // would be positive.
      
      if ( q > 0 )
        p = -p;
      else
        q = -q;

      // If x+p/q falls in [a,b] and is not too close to a and b, and isn't
      // too large, it is accepted.
      // If p/q is too large then the gold section procedure would
      // reduce [a,b] to larger extent.
      
      if (    (abs(p) < abs(new_step*q))
           && (p > q*(a-x+2*eps_act))
           && (p < q*(b-x-2*eps_act))  )
        new_step = p/q;
    }

    // Adjust the step to be not less than the tolerance.
    
    if ( abs(new_step) < eps_act )
      new_step =  new_step > 0 ? eps_act : -eps_act;

    // Obtain the next approximation to min
    // and reduce the encompassing interval.
    
    register const Real t = x + new_step;
    register const Real ft = f(t);

    if ( ft <= fx ) // t is a better approximation.
    {
      // Reduce the interval so that t would fall within it.
      
      ( (t < x) ? b : a ) = x;	

		// Assign the best approx to x.
      
      v = w;  w = x;  x = t;		
      fv=fw;  fw=fx;  fx=ft;
      
    }
    else // x remains the better approx.
    {
      // Reduce the interval encompassing x.
      
      ( (t < x) ? a : b ) = t;	
      
      if ( (ft <= fw) || (w==x) )
      {
        v = w;  w = t;
        fv=fw;  fw=ft;
      }
      else if ( (ft<=fv) || (v==x) || (v==w) )
      {
        v = t; fv = ft;
      } 
    }
  }

  cerr << "Warning: Brent minimiser did not reach requested accuracy "
       << "in interval [" << ax << "," << bx << "]. " << endl;

  return x;
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_minimum
//
//   Same as above, but can handle different intervals at the same time.
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_minimum(Function1D<Real>& f, vector<Real>& Ax,
                           vector<Real>& Bx, Real eps)
{
  if (Ax.size() != Bx.size())
  {
    cerr << "Error: Ax and Bx must be the same size." << endl;
    exit (-1);
  }

  vector<Real> minima;

  for (unsigned int i=0; i<Ax.size(); i++)
    minima.push_back(brent_minimum(f, Ax[i], Bx[i], eps));

  return minima;
}



/////////////////////////////////////////////////////////////////////////////
//
// bracket_all_minima
//
/////////////////////////////////////////////////////////////////////////////

void bracket_all_minima(Function1D<Real>& f, Real ax, Real bx,
                        vector<Real>& Ax,vector<Real>& Bx, Real dx,
                        int sec_level)
{
  Ax.clear(); Bx.clear();
  
  Real fx, previous_fx;
  Real fx_coarse, previous_fx_coarse;

  fx = previous_fx = fx_coarse = previous_fx_coarse = f(ax);

  bool decreasing = false, decreasing_coarse = false;

  long int iters = 0;
  int coarse_minima = 0;
  
  for (Real x=ax+dx; x<=bx; x+=dx)
  {
    fx = f(x);
    
    if ( decreasing && (fx > previous_fx) )
    {
      Ax.push_back(x-2*dx);
      Bx.push_back(x+2*dx);
    }
    
    decreasing = (fx < previous_fx);    
    previous_fx = fx;

    if (sec_level > 0)
    {
      if ( ((iters % int(pow(2.0,sec_level)))==0) || (x+dx>bx) )
      {
        fx_coarse = fx;
        
        if ( decreasing_coarse && (fx_coarse > previous_fx_coarse) )
          coarse_minima++;

        decreasing_coarse = (fx_coarse < previous_fx_coarse);
        previous_fx_coarse = fx_coarse;
      }
    }
  }

  if ( (sec_level > 0) && (int(Ax.size()) != coarse_minima) )
  {
    cout << "Warning: number of minima with coarse grid different than "
         << "with fine grid. " << endl;
    cout << "Step   dx: " << Ax.size() << endl;
    cout << "Step " << pow(2.0,sec_level) << "*dx: " << coarse_minima << endl;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// bracket_N_minima
//
/////////////////////////////////////////////////////////////////////////////

void bracket_N_minima(Function1D<Real>& f, Real ax, int N,
                      vector<Real>& Ax, vector<Real>& Bx,
                      Real dx, int sec_level)
{
  Ax.clear(); Bx.clear();

  Real fx, previous_fx;
  Real fx_coarse, previous_fx_coarse;

  fx = previous_fx = fx_coarse = previous_fx_coarse = f(ax);

  bool decreasing = false, decreasing_coarse = false;

  const long int MAXITER = 10000000L;
  long int iters         = 0;
       int fine_minima   = 0;
       int coarse_minima = 0;
  
  for (Real x=ax+dx; fine_minima<N; x+=dx, iters++)
  {
    fx = f(x);
    
    if ( decreasing && (fx > previous_fx) )
    {      
      Ax.push_back(x-2*dx);
      Bx.push_back(x+2*dx);
      fine_minima++;
    }

    decreasing = (fx < previous_fx);
    previous_fx = fx;

    if (sec_level > 0)
    {
      if (fine_minima==N)
        coarse_minima++; // Avoid problems at the end of the interval.
      else if ( (iters % int(pow(2.0,sec_level)))==0 )
      { 
        fx_coarse = fx;
        
        if ( decreasing_coarse && (fx_coarse > previous_fx_coarse) )
          coarse_minima++;
        
        decreasing_coarse = (fx_coarse < previous_fx_coarse);
        previous_fx_coarse = fx_coarse;
      }
    }  
  }

  if (iters >= MAXITER-1)
    cout << "Warning: maximum iterations " << MAXITER << " reached." << endl;

  if ( (sec_level > 0) && (int(Ax.size()) != coarse_minima) )
  {
    cout << "Warning: number of minima with coarse grid different than "
         << "with fine grid. " << endl;
    cout << "Step   dx: " << Ax.size() << endl;
    cout << "Step " << pow(2.0,sec_level) << "*dx: " << coarse_minima << endl;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_all_minima
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_all_minima(Function1D<Real>& f, Real ax, Real bx, Real dx,
                              Real eps, int sec_level)
{
  vector<Real> Ax, Bx;
  bracket_all_minima(f, ax, bx, Ax, Bx, dx, sec_level);
  return brent_minimum(f, Ax, Bx, eps);
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_N_minima
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_N_minima(Function1D<Real>& f, Real ax, int N, Real dx,
                            Real eps, int sec_level)
{
  vector<Real> Ax, Bx;
  bracket_N_minima(f, ax, N, Ax, Bx, dx, sec_level);
  return brent_minimum(f, Ax, Bx, eps); 
}



/////////////////////////////////////////////////////////////////////////////
//
// brent_refine_minima
//
/////////////////////////////////////////////////////////////////////////////

vector<Real> brent_refine_minima(Function1D<Real>& f, vector<Real>& x,
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
      = brent_all_minima(f, ax, bx, dx, eps, sec_level);

    for (unsigned int j=0; j<x_i.size(); j++)
      results.push_back(x_i[j]);
  }
  
  return results;
}



