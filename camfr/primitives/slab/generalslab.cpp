
/////////////////////////////////////////////////////////////////////////////
//
// File:     generalslab.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
//
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "generalslab.h"
#include "slabmatrixcache.h"
#include "isoslab/slabwall.h"
#include "../../math/calculus/quadrature/patterson_quad.h"

using std::vector;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// SlabGlobal
//
/////////////////////////////////////////////////////////////////////////////

SlabGlobal global_slab = {0.0, 0.0, NULL, NULL};



/////////////////////////////////////////////////////////////////////////////
//
// SlabFlux
//
/////////////////////////////////////////////////////////////////////////////

class SlabFlux : public RealFunction
{
  public:

    SlabFlux(const FieldExpansion& fe_) : fe(fe_) {}

    Real operator()(const Real& x)
    {
      counter++;
      Field f=fe.field(Coord(x,0,0));
      return real(f.E1*conj(f.H2) - f.E2*conj(f.H1));
    }

  protected:

    FieldExpansion fe;
};



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::~SlabImpl
//
/////////////////////////////////////////////////////////////////////////////

SlabImpl::~SlabImpl(){
  slabmatrix_cache.deregister(this);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::R_lower()
//
/////////////////////////////////////////////////////////////////////////////

Complex SlabImpl::R_lower() const
{
  SlabWall* l_wall = lowerwall ? lowerwall : global_slab.lowerwall;
  Complex R_lower = l_wall ? l_wall->get_R12() : -1.0;
  return R_lower;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::R_upper()
//
/////////////////////////////////////////////////////////////////////////////

Complex SlabImpl::R_upper() const
{
  SlabWall* u_wall = upperwall ? upperwall : global_slab.upperwall;
  Complex R_upper = u_wall ? u_wall->get_R12() : -1.0;
  return R_upper;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::S_flux()
//
/////////////////////////////////////////////////////////////////////////////

Real SlabImpl::S_flux(const FieldExpansion& f,
                      Real c1_start, Real c1_stop,
                      Real precision) const
{
  SlabFlux flux(f);
  return patterson_quad(flux, c1_start, c1_stop, precision);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::disc_intersect
//
//   Make sorted list of evaluation points for field cache.
//
/////////////////////////////////////////////////////////////////////////////

vector<Complex> SlabImpl::disc_intersect(const SlabImpl* medium_II) const
{
  vector<Complex> disc = discontinuities;

  disc.push_back(0.0);
  
  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-9);

  sort(disc.begin(), disc.end(), RealSorter());

  return disc;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::fill_field_cache
//
/////////////////////////////////////////////////////////////////////////////

// Dirty includes. Refactor this code after the design of moslab has settled.

#include "isoslab/slabmode.h"
#include "isoslab/slaboverlap.h"

void SlabImpl::fill_field_cache(SlabCache* cache, SlabImpl* medium_II,
                                const vector<Complex>& disc)
{  
  for (int i=1; i<=int(N()); i++)
  {
    const SlabMode* mode_I
      = dynamic_cast<const SlabMode*>(get_mode(i));
    
    const SlabMode* mode_II
      = dynamic_cast<const SlabMode*>(medium_II->get_mode(i));

    for (int k=0; k<int(disc.size()-1); k++)
    {
      const Coord lower(disc[k],  0,0, Plus);
      const Coord upper(disc[k+1],0,0, Min);

      Complex fw_I_l, bw_I_l, fw_II_l, bw_II_l;
      Complex fw_I_u, bw_I_u, fw_II_u, bw_II_u;

      mode_I ->forw_backw_at(lower, &fw_I_l,  &bw_I_l);    
      mode_I ->forw_backw_at(upper, &fw_I_u,  &bw_I_u);
      mode_II->forw_backw_at(lower, &fw_II_l, &bw_II_l);  
      mode_II->forw_backw_at(upper, &fw_II_u, &bw_II_u);

      cache->fw_l(1,i,k+1) = fw_I_l;
      cache->bw_l(1,i,k+1) = bw_I_l;
      cache->fw_u(1,i,k+1) = fw_I_u;
      cache->bw_u(1,i,k+1) = bw_I_u;

      cache->fw_l(2,i,k+1) = fw_II_l;
      cache->bw_l(2,i,k+1) = bw_II_l;
      cache->fw_u(2,i,k+1) = fw_II_u;
      cache->bw_u(2,i,k+1) = bw_II_u;
    }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

void SlabImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I, cMatrix* O_II_II)
{
  // Fill field cache.

  SlabImpl* medium_I  = this;
  SlabImpl* medium_II = dynamic_cast<SlabImpl*>(w);

  vector<Complex> disc = medium_I->disc_intersect(medium_II);

  const unsigned int N = global.N;
  SlabCache cache(N, disc.size()-1);
  medium_I->fill_field_cache(&cache, medium_II, disc);

  // Calculate overlap matrices (y-invariant case).
  
  if (global.polarisation != TE_TM)
  {
    for (int i=1; i<=int(N); i++)
      for (int j=1; j<=int(N); j++)
      {        
        (*O_I_II)(i,j) = overlap
          (dynamic_cast<const SlabMode*>(medium_I ->get_mode(i)),
           dynamic_cast<const SlabMode*>(medium_II->get_mode(j)),
           &cache, &disc, i, j, 1, 2);
      
        (*O_II_I)(i,j) = overlap
          (dynamic_cast<const SlabMode*>(medium_II->get_mode(i)),
           dynamic_cast<const SlabMode*>(medium_I ->get_mode(j)),
           &cache, &disc, i, j, 2, 1);

        if (O_I_I) (*O_I_I)(i,j) = overlap
          (dynamic_cast<const SlabMode*>(medium_I ->get_mode(i)),
           dynamic_cast<const SlabMode*>(medium_I ->get_mode(j)),
           &cache, &disc, i, j, 1, 1);

        if (O_II_II) (*O_II_II)(i,j) = overlap
          (dynamic_cast<const SlabMode*>(medium_II->get_mode(i)),
           dynamic_cast<const SlabMode*>(medium_II->get_mode(j)),
           &cache, &disc, i, j, 2, 2);
      }

    return;
  }

  // Calculate overlap matrices for off-angle incidence.

  if (!O_I_I)
  { 
    py_error(
     "Internal error: non-orthogonality of modes not taken into account.");
    exit (-1);
  }

  cVector sin_I (N,fortranArray), cos_I (N,fortranArray);
  cVector sin_II(N,fortranArray), cos_II(N,fortranArray);

  for (int i=1; i<=int(N); i++)
  {
    SlabMode* mode_I  = dynamic_cast<SlabMode*>(medium_I ->get_mode(i));
    SlabMode* mode_II = dynamic_cast<SlabMode*>(medium_II->get_mode(i));   
    
    sin_I (i) = mode_I ->get_sin(); cos_I (i) = mode_I ->get_cos();
    sin_II(i) = mode_II->get_sin(); cos_II(i) = mode_II->get_cos();

    //std::cout << i <<cos_I (i) << cos_I(i)*mode_I->get_kz0() 
    //          << cos_II (i) << cos_II(i)*mode_II->get_kz0() 
    //          << sqrt(cos_I(i)) << sqrt(cos_II(i)) << std::endl;
  }

  // Calculate O_I_II and O_II_I.
  
  const int n = int(N/2);

  const OverlapMatrices* m 
    = slabmatrix_cache.get_matrices(medium_I, medium_II, &cache, &disc);

  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
    {
      (*O_I_II)(i,    j) =   m->TE_TE      (1,i,j) * cos_I (  i);
      (*O_I_II)(i,  n+j) =   0.0;
      (*O_I_II)(n+i,  j) =   m->Ex_Hz_cross(1,i,j) * sin_II(  j)
                           - m->Ez_Hx_cross(1,i,j) * sin_I (n+i);   
      (*O_I_II)(n+i,n+j) =   m->TM_TM      (1,i,j) * cos_II(n+j);

      (*O_II_I)(i,    j) =   m->TE_TE      (2,i,j) * cos_II(  i);
      (*O_II_I)(i,  n+j) =   0.0;
      (*O_II_I)(n+i,  j) =   m->Ex_Hz_cross(2,i,j) * sin_I (  j)
                           - m->Ez_Hx_cross(2,i,j) * sin_II(n+i);
      (*O_II_I)(n+i,n+j) =   m->TM_TM      (2,i,j) * cos_I (n+j);
    }

  for (int i=1; i<=N; i++)
    for (int j=1; j<=N; j++)
    {
      (*O_I_II)(i,j) /= sqrt(cos_I (i)) * sqrt(cos_II(j));
      (*O_II_I)(i,j) /= sqrt(cos_II(i)) * sqrt(cos_I (j));
    }
  

  // Calculate O_I_I and O_II_II.

  *O_I_I = 0.0; *O_II_II = 0.0;
  
  for (int i=1; i<=N; i++)
    (*O_I_I)(i,i) = (*O_II_II)(i,i) = 1.0;  

  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
    {
      (*O_I_I)(n+i,j)   = (  m->Ex_Hz_self(1,i,j) * sin_I(  j)
                           - m->Ez_Hx_self(1,i,j) * sin_I(n+i))
                             / sqrt(cos_I(n+i)) / sqrt(cos_I(j));
      
      (*O_II_II)(n+i,j) = (  m->Ex_Hz_self(2,i,j) * sin_II(  j)
                           - m->Ez_Hx_self(2,i,j) * sin_II(n+i))
                             / sqrt(cos_II(n+i)) / sqrt(cos_II(j));
    }
}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: OverlapFunction
//
//  Integrandum for overlap integral of a mode with an arbitrary field 
//  profile described by a function f. For TE, f is E2, for TM f is H2.
//
/////////////////////////////////////////////////////////////////////////////

class OverlapFunction : public ComplexFunction
{
  public:

    OverlapFunction(SlabMode* m_, ComplexFunction* f_) : m(m_), f(f_) {}

    Complex operator()(const Complex& x)
    {
      counter++;

      if (m->pol == TE)
        return -(*f)(x) * m->field(Coord(x,0,0)).H1;
      else
        return  (*f)(x) * m->field(Coord(x,0,0)).E1;
    }

  protected:

    SlabMode* m;
    ComplexFunction* f;
};



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::expand_field
//
//  Returns the expansion coefficients of a general field profile f.
//  Eps is the precision for the numerical integration.
// 
/////////////////////////////////////////////////////////////////////////////

cVector SlabImpl::expand_field(ComplexFunction* f, Real eps)
{
  find_modes();

  cVector coef(global.N, fortranArray); 
  coef = 0.0;
  
  vector<Complex> disc = discontinuities;
  disc.insert(disc.begin(), 0.0);

  for(int m=1; m<=global.N; m++)
  {
    SlabMode* mode = dynamic_cast<SlabMode*>(get_mode(m));

    OverlapFunction o(mode, f);
    
    Wrap_real_to_real f_re(o);
    Wrap_real_to_imag f_im(o);

    // Speed up convergence by splitting the integrals.

    for(int k=0; k<int(disc.size()-1); k++)
    {     
      Real begin = real(disc[k]);
      Real end   = real(disc[k+1]);
      coef(m) +=   patterson_quad(f_re, begin, end, eps);
      coef(m) += I*patterson_quad(f_im, begin, end, eps);
    }
  }

  return coef;
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab::Slab
//
/////////////////////////////////////////////////////////////////////////////

Slab::Slab(const Expression& ex, int M_series)
{
  Expression e = ex.flatten();
  
  if (e.get_size() == 1)
  {
    Material* m = dynamic_cast<Material*>(e.get_term(0)->get_mat());
    Complex   d = e.get_term(0)->get_d() 
      + I*global_slab.lower_PML + I*global_slab.upper_PML;
    
    if (!m)
    {
      py_error("Error: expression contains non-material term.");
      return;
    }
    
    s = new UniformSlab(d, *m);
  }
  else
    s = new Slab_M(e,M_series);

  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab::Slab
//
/////////////////////////////////////////////////////////////////////////////

Slab::Slab(const Term& t, int M_series)
{
  Material* m = dynamic_cast<Material*>(t.get_mat());
  Complex   d = t.get_d() + I*global_slab.lower_PML + I*global_slab.upper_PML;

  if (!m)
  {
    py_error("Error: expression contains non-material term.");
    return;
  }

  s = new UniformSlab(d, *m);  

  uniform = s->is_uniform();
  core = s->get_core();
}
