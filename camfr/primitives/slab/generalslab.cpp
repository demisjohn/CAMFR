
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

// Work around MS linker bug.

#ifdef _WIN32
#include "slabmatrixcache.cpp"
#else
#include "slabmatrixcache.h"
#endif

#include "generalslab.h"
#include "../../math/calculus/quadrature/patterson_quad.h"

using std::vector;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// SlabGlobal
//
/////////////////////////////////////////////////////////////////////////////

SlabGlobal global_slab = {NULL, NULL, 0.0};



/////////////////////////////////////////////////////////////////////////////
//
// SlabFlux
//
/////////////////////////////////////////////////////////////////////////////

class SlabFlux : public RealFunction
{
  public:

    SlabFlux(const FieldExpansion& fe_) : fe(fe_) 
    {
      PML = dynamic_cast<SlabImpl*>(fe.wg)->get_imag_start_thickness();
    }

    Real operator()(const Real& x)
    {
      counter++;
      Field f=fe.field(Coord(x+PML*I,0,0));
      return real(f.E1*conj(f.H2) - f.E2*conj(f.H1));
    }

  protected:

    FieldExpansion fe;
    Real PML;
};



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::~SlabImpl
//
/////////////////////////////////////////////////////////////////////////////

SlabImpl::~SlabImpl() 
{
  slabmatrix_cache.deregister(this);
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
// signedsqrt2
//
//   Square root with branch cut at 45 degrees.
//
/////////////////////////////////////////////////////////////////////////////

Complex signedsqrt2(const Complex& kz2)
{
  Complex new_kz = sqrt(kz2);

  if (imag(new_kz) > 0)
    new_kz = -new_kz;

  if (abs(imag(new_kz)) < abs(real(new_kz)))
    if (real(new_kz) < 0)
      new_kz = -new_kz;

  return new_kz;
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

// Dirty includes. Refactor this code after the design of moslab has settled.

#include "isoslab/slabmode.cpp"

// Work around MS linker bug.

#ifdef _WIN32
#include "isoslab/slaboverlap.cpp"
#else
#include "isoslab/slaboverlap.h"
#endif

void SlabImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I, cMatrix* O_II_II)
{ 
  SlabImpl* medium_I  = this;
  SlabImpl* medium_II = dynamic_cast<SlabImpl*>(w);

  // Make sorted list of evaluation points for field cache.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);
  
  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  sort(disc.begin(), disc.end(), RealSorter());

  // Fill field cache.

  const unsigned int N = global.N;

  SlabCache cache(N, disc.size()-1);
  
  for (int i=1; i<=int(N); i++)
  {
    const SlabMode* mode_I
      = dynamic_cast<const SlabMode*>(medium_I ->get_mode(i));
    
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

      cache.fw_l(1,i,k+1) = fw_I_l;
      cache.bw_l(1,i,k+1) = bw_I_l;
      cache.fw_u(1,i,k+1) = fw_I_u;
      cache.bw_u(1,i,k+1) = bw_I_u;

      cache.fw_l(2,i,k+1) = fw_II_l;
      cache.bw_l(2,i,k+1) = bw_II_l;
      cache.fw_u(2,i,k+1) = fw_II_u;
      cache.bw_u(2,i,k+1) = bw_II_u;
    }
  }

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
     "Internal error: non-orthogonality of modes not taken into account");
    exit (-1);
  }

  cVector sin_I (N,fortranArray), cos_I (N,fortranArray);
  cVector sin_II(N,fortranArray), cos_II(N,fortranArray);

  for (int i=1; i<=int(N); i++)
  {
    Complex kz0_I  
      = dynamic_cast<SlabMode*>(medium_I ->get_mode(i))->get_kz0();
    Complex kz0_II 
      = dynamic_cast<SlabMode*>(medium_II->get_mode(i))->get_kz0();
    
    sin_I (i) = global_slab.beta / kz0_I;
    sin_II(i) = global_slab.beta / kz0_II;

    // Note that cos needs to be in sync with the ones in slabmode.cpp.

    cos_I (i) = signedsqrt2(1.0 - sin_I (i) * sin_I (i));
    cos_II(i) = signedsqrt2(1.0 - sin_II(i) * sin_II(i));
  }

  // Calculate O_I_II and O_II_I.

  const int n = int(N/2);

  const OverlapMatrices* m 
    = slabmatrix_cache.get_matrices(medium_I, medium_II, &cache, &disc);

  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
    {
      (*O_I_II)(i,    j) =   m->TE_TE_I_II(i,j) * cos_I (  i);
      (*O_I_II)(i,  n+j) =   0.0;
      (*O_I_II)(n+i,  j) =   m->Ex_Hz_I_II(i,j) * sin_II(  j)
                           - m->Ez_Hx_I_II(i,j) * sin_I (n+i);   
      (*O_I_II)(n+i,n+j) =   m->TM_TM_I_II(i,j) * cos_II(n+j);

      (*O_II_I)(i,    j) =   m->TE_TE_II_I(i,j) * cos_II(  i);
      (*O_II_I)(i,  n+j) =   0.0;
      (*O_II_I)(n+i,  j) =   m->Ex_Hz_II_I(i,j) * sin_I (  j)
                           - m->Ez_Hx_II_I(i,j) * sin_II(n+i);
      (*O_II_I)(n+i,n+j) =   m->TM_TM_II_I(i,j) * cos_I (n+j);  
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
      (*O_I_I)(n+i,j)   = (  m->Ex_Hz_I_I(i,j) * sin_I(  j)
                           - m->Ez_Hx_I_I(i,j) * sin_I(n+i)) 
                              / sqrt(cos_I(n+i)) / sqrt(cos_I(j));

      (*O_II_II)(n+i,j) = (  m->Ex_Hz_II_II(i,j) * sin_II(  j)
                           - m->Ez_Hx_II_II(i,j) * sin_II(n+i))
                              / sqrt(cos_II(n+i)) / sqrt(cos_II(j));
    }
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab::Slab
//
/////////////////////////////////////////////////////////////////////////////

Slab::Slab(const Expression& ex)
{
  Expression e = ex.flatten();
  
  if (e.get_size() == 1)
  {
    Material* m = dynamic_cast<Material*>(e.get_term(0)->get_mat());
    Complex   d = e.get_term(0)->get_d();
    
    if (!m)
    {
      py_error("Error: expression contains non-material term.");
      return;
    }
    
    s = new UniformSlab(d, *m);
  }
  else
    s = new Slab_M(e);

  uniform = s->is_uniform();
  core = s->get_core();
}



/////////////////////////////////////////////////////////////////////////////
//
// Slab::Slab
//
/////////////////////////////////////////////////////////////////////////////

Slab::Slab(const Term& t)
{
  Material* m = dynamic_cast<Material*>(t.get_mat());
  Complex   d = t.get_d();

  if (!m)
  {
    py_error("Error: expression contains non-material term.");
    return;
  }

  s = new UniformSlab(d, *m);

  uniform = s->is_uniform();
  core = s->get_core();
}
