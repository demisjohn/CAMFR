
/////////////////////////////////////////////////////////////////////////////
//
// File:     contour.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010322
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "../../../util/cvector.h"
#include "patterson_z_n.h"
#include "contour.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// Description of the paths of the main contour and the subcontours.
// For sub_segments and sub_signs, the first index is a Subcontour index.
//
// Points on the contour are labelled as follows:
//
//                  tl  tc  tr
//                  cl  cc  cr
//                  bl  bc  br
//
/////////////////////////////////////////////////////////////////////////////

Segment Contour::main_segments[8]
  = {br_cr, cr_tr, tr_tc, tc_tl, tl_cl, cl_bl, bl_bc, bc_br};
int     Contour::main_signs[8]
  = {+1,    +1,    +1,    +1,    +1,    +1,    +1,    +1};

Segment Contour::sub_segments[4][4] = { {cc_tc, tc_tl, tl_cl, cc_cl},
                                        {cr_tr, tr_tc, cc_tc, cc_cr},
                                        {cc_bc, cc_cl, cl_bl, bl_bc},
                                        {br_cr, cc_cr, cc_bc, bc_br} };

int     Contour::sub_signs[4][4]    = { {+1,    +1,    +1,    -1},
                                        {+1,    +1,    -1,    +1},
                                        {-1,    +1,    +1,    +1},
                                        {+1,    -1,    +1,    +1} };



/////////////////////////////////////////////////////////////////////////////
//
// Contour::Contour
//
/////////////////////////////////////////////////////////////////////////////

Contour::Contour(const Complex& bottom_left, const Complex& top_right,
                 ComplexFunction& f_, unsigned int M_,
                 Real eps_=1e-4, Real mu_=1e-4, unsigned int max_k_=8)
  : bl(bottom_left), tr(top_right),
    f(&f_), M(M_), eps(eps_), mu(mu_), max_k(max_k_)
{
  for (unsigned int i=0; i<12; i++)
    know_integrals[i] = false;

  create_internal_points();
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::create_internal_points
//
/////////////////////////////////////////////////////////////////////////////

void Contour::create_internal_points()
{  
  // Construct egde points on contours.
  
  tl = real(bl) + I*imag(tr);
  br = real(tr) + I*imag(bl);

  tc = (tl + tr) / 2.;
  bc = (bl + br) / 2.;
  cl = (bl + tl) / 2.;
  cr = (br + tr) / 2.;

  cc = (cl + cr) / 2.;
  
  // Create lookup table

  begin[br_cr] = &br; end[br_cr] = &cr;
  begin[cr_tr] = &cr; end[cr_tr] = &tr;
  begin[tr_tc] = &tr; end[tr_tc] = &tc;
  begin[tc_tl] = &tc; end[tc_tl] = &tl;  
  begin[tl_cl] = &tl; end[tl_cl] = &cl;
  begin[cl_bl] = &cl; end[cl_bl] = &bl;
  begin[bl_bc] = &bl; end[bl_bc] = &bc;
  begin[bc_br] = &bc; end[bc_br] = &br;
  begin[cc_tc] = &cc; end[cc_tc] = &tc;
  begin[cc_cr] = &cc; end[cc_cr] = &cr;
  begin[cc_bc] = &cc; end[cc_bc] = &bc;
  begin[cc_cl] = &cc; end[cc_cl] = &cl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::Contour
//
/////////////////////////////////////////////////////////////////////////////

Contour::Contour(const Contour& c)
{
  copy_from(c);
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::operator=
//
/////////////////////////////////////////////////////////////////////////////
    
Contour& Contour::operator=(const Contour& c)
{
  if (this == &c)
    return *this;

  copy_from(c);

  return *this;
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::copy_from
//
/////////////////////////////////////////////////////////////////////////////

void Contour::copy_from(const Contour& c)
{
  bl    = c.bl; 
  tr    = c.tr;
  f     = c.f;
  M     = c.M;
  eps   = c.eps;
  mu    = c.mu;
  max_k = c.max_k;
  
  for (unsigned int i=0; i<12; i++)
  {
    know_integrals[i] = c.know_integrals[i];
         integrals[i] = c.     integrals[i];
  }

  create_internal_points();
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::encloses
//
/////////////////////////////////////////////////////////////////////////////

bool Contour::encloses(const Complex& z) const
{  
  return ( (real(bl)<=real(z)) && (real(z)<=real(tr)) &&
           (imag(bl)<=imag(z)) && (imag(z)<=imag(tr)) );
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::subcontour
//
/////////////////////////////////////////////////////////////////////////////

const Contour Contour::subcontour(Subcontour s) const
{
  switch(s)
  {
    case top_left:
      return Contour(cl, tc, *f, M, eps, mu, max_k);

    case top_right:
      return Contour(cc, tr, *f, M, eps, mu, max_k);

    case bottom_left:
      return Contour(bl, cc, *f, M, eps, mu, max_k);

    case bottom_right:
      return Contour(bc, cr, *f, M, eps, mu, max_k);
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::get_integrals
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Contour::get_integrals(Segment segment) const
{
  if (know_integrals[segment] == true)
    return integrals[segment];

  integrals[segment]
    = patterson_quad_z_n(*f,*begin[segment],*end[segment],M,eps,mu,max_k);
  know_integrals[segment] = true;
  
  return integrals[segment];
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::path_integrals
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> Contour::path_integrals
  (unsigned int no_segments, Segment segments[], int signs[]) const
{
  // Check the smallest index of M (this is an adaptive parameter that can
  // be decreased in the patterson quadrature.

  unsigned int M_min = get_integrals(segments[0]).size();
  for (unsigned int i=1; i<no_segments; i++)
    if (get_integrals(segments[i]).size() < M_min)
      M_min = get_integrals(segments[i]).size();

  // Accumulate result.

  vector<Complex> result;
  for (unsigned int m=0; m<M_min; m++)
  {
    Complex result_m = 0.0;
    for (unsigned int i=0; i<no_segments; i++)
      result_m += Real(signs[i])*get_integrals(segments[i])[m];
    
    result.push_back(result_m);
  }
  
  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::adjacent_r
//
/////////////////////////////////////////////////////////////////////////////

Contour Contour::adjacent_r() const
{
  Contour r(br, tr + (tr-tl), *f, M, eps, mu, max_k);

  r.set_integrals(tl_cl, -get_integrals(cr_tr));
  r.set_integrals(cl_bl, -get_integrals(br_cr));

  cout << "Extending " << bl << tr << endl;
    cout << "to " << r.get_bottom_left()
         << r.get_top_right() << endl;

  return r;
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::adjacent_ur
//
/////////////////////////////////////////////////////////////////////////////

vector<Contour> Contour::adjacent_ur() const
{
  // Create adjacent contours.

  Contour  r(br, tr + (tr-tl), *f, M, eps, mu, max_k);
  Contour  u(tl, tr + (tr-br), *f, M, eps, mu, max_k);
  Contour ur(tr, tr + (tr-bl), *f, M, eps, mu, max_k);

  // Precalculate shared line segments.

   r.set_integrals(tl_cl,   -get_integrals(cr_tr));
   r.set_integrals(cl_bl,   -get_integrals(br_cr));

  ur.set_integrals(tl_cl, -u.get_integrals(cr_tr));
  ur.set_integrals(cl_bl, -u.get_integrals(br_cr));

   u.set_integrals(bl_bc,   -get_integrals(tc_tl));
   u.set_integrals(bc_br,   -get_integrals(tr_tc));

  ur.set_integrals(bl_bc, -r.get_integrals(tc_tl));
  ur.set_integrals(bc_br, -r.get_integrals(tr_tc));
  
  // Return adjacent contours.

  vector<Contour> new_contours;

  new_contours.push_back(ur);
  new_contours.push_back(u);
  new_contours.push_back(r);

  cout << "Extending " << bl << tr << endl;
  for (unsigned int i=0; i<new_contours.size(); i++)
    cout << "to " << i << new_contours[i].get_bottom_left()
         << new_contours[i].get_top_right() << endl;

  return new_contours;
}



/////////////////////////////////////////////////////////////////////////////
//
// Contour::double_ur
//
/////////////////////////////////////////////////////////////////////////////

Contour Contour::double_ur() const
{
  Contour ur(bl, tr + (tr-bl), *f, M, eps, mu, max_k);

  ur.set_integrals(cc_bc, -get_integrals(br_cr) - get_integrals(cr_tr));
  ur.set_integrals(cc_cl,  get_integrals(tr_tc) + get_integrals(tc_tl));
  ur.set_integrals(cl_bl,  get_integrals(tl_cl) + get_integrals(cl_bl));
  ur.set_integrals(bl_bc,  get_integrals(bl_bc) + get_integrals(bc_br));

  cout << "Enlarging initial contour to" << ur.get_bottom_left()
       << ur.get_top_right() << endl;

  return ur;
}
