
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
#include "../../math/calculus/quadrature/patterson_quad.h"
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

    SlabFlux(const FieldExpansion& fe_) : fe(fe_) {}

    Real operator()(const Real& x)
    {
      counter++;
      Field f=fe.field(Coord(x,0,0));
      return (global.polarisation == TE) ? real(-f.E2*conj(f.H1))
                                         : real( f.E1*conj(f.H2));
    }

  protected:

    FieldExpansion fe;
};



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::S_flux()
//
/////////////////////////////////////////////////////////////////////////////

Real SlabImpl::S_flux(const FieldExpansion& f,
                      Real c1_start, Real c1_stop,
                      Real precision = 1e-10) const
{
  SlabFlux flux(f);
  return patterson_quad(flux, c1_start, c1_stop, precision);
}



/////////////////////////////////////////////////////////////////////////////
//
// SlabImpl::calc_overlap_matrices()
//
/////////////////////////////////////////////////////////////////////////////

// Dirty includes. Refactor this code after the design of moslab has settled.

#include "isoslab/slaboverlap.h"

void SlabImpl::calc_overlap_matrices
  (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
   cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL)
{  
  const SlabImpl* medium_I  = this;
  const SlabImpl* medium_II = dynamic_cast<const SlabImpl*>(w);

  // Make sorted list of evaluation points for field cache.

  vector<Complex> disc = medium_I->discontinuities;

  disc.push_back(0.0);
  
  for (unsigned int k=0; k<medium_II->discontinuities.size(); k++)
    disc.push_back(medium_II->discontinuities[k]);

  remove_copies(&disc, 1e-6);

  sort(disc.begin(), disc.end(), RealSorter());

  // Fill field cache.

  SlabCache cache(global.N, disc.size()-1);
  
  for (int i=1; i<=int(global.N); i++)
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

  // Calculate overlap matrices.  

  for (int i=1; i<=int(global.N); i++)
    for (int j=1; j<=int(global.N); j++)
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
      cerr << "Error: expression contains non-material term." << endl;
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
    cerr << "Error: expression contains non-material term." << endl;
    return;
  }

  s = new UniformSlab(d, *m);

  uniform = s->is_uniform();
  core = s->get_core();
}
