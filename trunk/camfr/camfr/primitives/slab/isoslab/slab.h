
/////////////////////////////////////////////////////////////////////////////
//
// File:     slab.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
//
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SLAB_H
#define SLAB_H

#include "string"
#include "algorithm"
#include "../generalslab.h"
#include "../../../waveguide.h"
#include "../../../expression.h"
#include "../../../util/index.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Slab_M
//
//   Stack of M uniform infinite planar layers sandwiched between two metal
//   plates.
//   The layer interfaces are parallel the the z-axis, the modes propagate
//   in the z-direction, the left plate is at x=0 and the structure is
//   invariant in the y-direction.
//  
/////////////////////////////////////////////////////////////////////////////

class Slab_M : public SlabImpl
{ 
  public:

    Slab_M() {}
    Slab_M(const Expression& expression);  
    Slab_M(const Slab_M& slab);
    ~Slab_M() {}

    Slab_M& operator= (const Slab_M& slab);
    bool  operator==(const Waveguide& w) const;

    vector<Material*> get_materials() const
      {return materials;}
    
    bool contains(const Material& m) const
      {return find(materials.begin(), materials.end(), &m) != materials.end();}

    bool no_gain_present() const;

    Complex eps_at(const Coord& coord) const
      {return materials
         [index_lookup(coord.c1, coord.c1_limit, discontinuities)]->eps();}
    
    Complex mu_at(const Coord& coord) const
      {return materials
         [index_lookup(coord.c1, coord.c1_limit, discontinuities)]->mu();}

    Complex get_width() const
      {return discontinuities[discontinuities.size()-1];}
    
    void find_modes();

  protected:

    void find_modes_from_scratch_by_ADR();
    void find_modes_from_scratch_by_track();
    void find_modes_by_sweep();

    vector<Complex> params; // Last parameters of dispersion relation.
    
    vector<Material*> materials;
    vector<Complex>   thicknesses;

    friend class Slab_M_Mode;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: UniformSlab
//  
/////////////////////////////////////////////////////////////////////////////

class UniformSlab : public SlabImpl
{ 
  public:

    UniformSlab(const Complex& width, Material& core_)
      {uniform = true; core = &core_; discontinuities.push_back(width);}

    UniformSlab(const Waveguide_length& wg_l)
      {core = wg_l.wg->get_core(); discontinuities.push_back(wg_l.d);}

    ~UniformSlab() {}
    
    bool operator==(const Waveguide& w) const
      {return *core == *(w.get_core());}
    
    vector<Material*> get_materials() const
      {vector<Material*> m; m.push_back(core); return m;}
    
    bool contains(const Material& m) const
      {return *core == m;}

    bool no_gain_present() const
      {return core->no_gain_present();}

    Complex eps_at(const Coord& coord) const
      {return core->eps();}
        
    Complex mu_at(const Coord& coord) const
      {return core->mu();}

    Complex get_width() const
      {return discontinuities[0];}

    void find_modes();

  protected:

    void find_modes_single_pol();
};



#endif
