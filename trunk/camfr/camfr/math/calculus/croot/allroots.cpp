
/////////////////////////////////////////////////////////////////////////////
//
// File:     allroots.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010322
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "../../../util/cvector.h"
#include "../../linalg/linalg.h"
#include "../polyroot/polyroot.h"
#include "allroots.h"
#include "mueller.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

/////////////////////////////////////////////////////////////////////////////
//
// Maximum number of roots in a contour before we invoke subdivision.
//
/////////////////////////////////////////////////////////////////////////////

const unsigned int N_max = 10;



/////////////////////////////////////////////////////////////////////////////
//
// roots_contour
//
//   Low level routine, used by allroots.
//
//   Locate at most N roots of the function f in a rectangular contour.
//   starting from the integrals z^n/f(z). The value of the integrals are
//   passed separately, rather than by calling contour.contour_integrals,
//   to allow reuse of calculations from mother contour.
//
//   Does not used adaptive subdivision, such that there is no guarantee
//   that all the roots are found.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> roots_contour(const Contour& contour,
                              const vector<Complex>& integrals)
{
  const Complex bl = contour.get_bottom_left();
  const Complex tr = contour.get_top_right();
  
  // Calculate G coefficients and adjust N because round-off detection
  // might have changed the number of coefficients.

  vector<Complex> G = integrals / 2. / pi / I;

  unsigned int N = (unsigned int)(ceil(G.size() / 2.0));

  // Check if all G's are zero.

  bool all_zeros = true;
  for (unsigned int i=0; i<=2*N-1; i++)
    if (abs(G[i]) > 1e-10)
      all_zeros = false;

  if (all_zeros)
  {
    vector<Complex> empty;
    return empty;
  }

  // Construct linear system.

  cMatrix A(N,N,fortranArray), B(N,1,fortranArray), c(N,1,fortranArray);

  for (int i=1; i<=N; i++)
    for (int j=1; j<=N; j++)
      A(i,j) = G[i+j-2];

  for (int i=1; i<=N; i++)
    B(i,1) = -G[N+i-1];

  // Solve linear system and exploit the symmetry. Note: the fact that this
  // is a Hankel matrix is not used, since speedups would probably be minor
  // due to the small dimension.

  c = solve_sym(A,B);

  // Find roots of polynomial.

  vector<Complex> coeff;
  coeff.push_back(1.0);
  for (int i=N; i>=1; i--)
    coeff.push_back(c(i,1));

  vector<Complex> roots1 = polyroot(coeff);

  // Refine them using mueller solver with deflation.

  vector<Complex> roots2, deflate;
  for (unsigned int i=0; i<roots1.size(); i++)
  {
    // If the estimate is far outside the contour, don't bother.

    if ( abs(roots1[i] - (bl+tr)/2.) > 1000*abs(bl-tr) )
      continue;
    
    bool error = false;
    
    Complex root = mueller(*contour.get_f(), roots1[i], roots1[i]+.001,
                           1e-14, &deflate, 50, &error);

    if (error)
      continue;

    // If the estimate converges to a root outside of the contour,
    // deflate it and try again.

    if (!contour.encloses(root))
    {      
      vector<Complex> outside_roots = mueller_multiple(*contour.get_f(), root);
      deflate.insert(deflate.end(),outside_roots.begin(),outside_roots.end());
      
      root = mueller(*contour.get_f(), roots1[i], roots1[i]+.001,
                     1e-14, &deflate, 50, &error);
    }
    
    if (contour.encloses(root))
    {      
      vector<Complex> multiple_roots 
        = mueller_multiple(*contour.get_f(),root,1e-14,&deflate);

      for (unsigned int j=0; j<multiple_roots.size(); j++)
        if (contour.encloses(multiple_roots[j]))
        {
          roots2.push_back(multiple_roots[j]);
          deflate.push_back(multiple_roots[j]);
        }
    }
    
  }
  
  return roots2;      
}



/////////////////////////////////////////////////////////////////////////////
//
// allroots_contour
//
//   Low level routine, used by allroots.
//
//   Recursively subdivides a contour to locate all roots of a function f
//   inside this contour.
//
/////////////////////////////////////////////////////////////////////////////

struct Sorter
{
    bool operator()(const Complex& a, const Complex& b)
       {return (real(a) < real(b));}
};

vector<Complex> allroots_contour(const Contour& c) 
{ 
  // Find roots in contour and subcontours.
  
  vector<Complex> roots = roots_contour(c, c.contour_integrals());

  bool OK = true;

  vector<Complex> subroots;
  for (unsigned int i_=0; i_<4; i_++)
  {
    Subcontour i = Subcontour(i_);
    vector<Complex> subroots_i
      = roots_contour(c.subcontour(i), c.subcontour_integrals(i));

    if (subroots_i.size() > 0.8*N_max)
      OK = false;
    
    subroots.insert(subroots.end(), subroots_i.begin(), subroots_i.end());
  }

  // If roots contains the same values as subroots, we are done.

  if (OK && (roots.size() == subroots.size()) )
  {
    Sorter sorter;

    sort(   roots.begin(),    roots.end(), sorter);
    sort(subroots.begin(), subroots.end(), sorter);

    bool equal = true;
    for (unsigned int i=0; i<roots.size(); i++)   
      if (abs(roots[i]-subroots[i]) > 1e-10)
        equal = false;
    
    if (equal)
      return roots;
  }

  // Recatch case of single lost zero. Heuristic to avoid having to
  // increase the precision unduly.

  if ((roots.size() == 1) && (subroots.size() == 0))
    return roots;
  
  // Else, recursively subdivide contour further.
  
  roots.clear();
  for (unsigned int i_=0; i_<4; i_++)
  {
    Subcontour i = Subcontour(i_);
    vector<Complex> subroots_i = allroots_contour(c.subcontour(i));
    roots.insert(roots.end(), subroots_i.begin(), subroots_i.end());
  }

  return roots;
}



/////////////////////////////////////////////////////////////////////////////
//
// allroots
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> allroots
  (ComplexFunction& f, const Complex& bottom_left, const Complex& top_right,
   Real eps, Real mu, unsigned int max_k)
{ 
  Contour contour(bottom_left, top_right, f, 2*N_max-1, eps, mu, max_k);
  return allroots_contour(contour);
}



/////////////////////////////////////////////////////////////////////////////
//
// N_roots
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> N_roots(ComplexFunction& f, unsigned int N,
                        const Complex& bottom_left, const Complex& top_right,
                        Real eps, Real mu, unsigned int max_k, 
                        ExpandDirection dir)
{
  // Enlarge initial contour if needed.

  Contour c0(bottom_left, top_right, f, 2*N-1, eps, mu, max_k);

  vector<Complex> roots
    = allroots(f, c0.get_bottom_left(), c0.get_top_right(), eps, mu, max_k);
  
  while (roots.size() == 0)
  {    
    Contour c = c0.double_ur();
    roots 
      = allroots(f, c.get_bottom_left(), c.get_top_right(), eps, mu, max_k);
   
    c0 = c;
  }

  if (roots.size() >= N)
    return roots;

  // Set up the contour stack with new contours.

  ExpandDirection true_dir = r;
  for (unsigned int i=0; i<roots.size(); i++)
  {
    if ( (dir == ur) && (imag(roots[i]) > imag(c0.get_center())) )
      true_dir = ur;
    else if ( (dir == dr) && (imag(roots[i]) < imag(c0.get_center())) )
      true_dir = dr;
  }

  vector<Contour> contour_stack;
  if (true_dir == r)
    contour_stack.push_back(c0.adjacent_r());
  else
  {
    vector<Contour> new_c = (true_dir == ur) ? c0.adjacent_ur() 
                                             : c0.adjacent_dr();

    for (unsigned int i=0; i<new_c.size(); i++)
      contour_stack.push_back(new_c[i]);
  }

  // Take additional contours into account until we have sufficient modes.

  int iters = 0;
  while (!contour_stack.empty())
  { 
    iters++;
    
    Contour c = contour_stack.front();
    contour_stack.erase(contour_stack.begin());
    
    vector<Complex> roots_c
      = allroots(f, c.get_bottom_left(), c.get_top_right(), eps, mu, max_k);

    ExpandDirection true_dir = r;
    for (unsigned int i=0; i<roots_c.size(); i++)
    {
      if ( (dir == ur) && (imag(roots_c[i]) > imag(c0.get_center())) )
        true_dir = ur;
      else if ( (dir == dr) && (imag(roots_c[i]) < imag(c0.get_center())) )
        true_dir = dr;
    }

    if (roots_c.size() == 0)
      true_dir = dir;

    roots.insert(roots.end(), roots_c.begin(), roots_c.end());

    if (    (roots.size() < N)
         && (!roots_c.empty() || (roots_c.empty() && contour_stack.empty())) )
    {

      if (true_dir == r)
        contour_stack.push_back(c.adjacent_r());
      else
      {
        vector<Contour> new_c = (true_dir == ur) ? c.adjacent_ur() 
                                                 : c.adjacent_dr();

        for (unsigned int i=0; i<new_c.size(); i++)
          contour_stack.push_back(new_c[i]);
      }
    }

    if (iters == 200)
    {
      py_print("Warning: maximum number of iterations reached.");
      py_print(
       "Try increasing real and/or imaginary part of set_C_upperright.");
      
      return roots;
    }
    
  }
  
  return roots;
}
