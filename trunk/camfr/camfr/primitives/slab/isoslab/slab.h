
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
    Slab_M(const Expression& expression, int M_series=0);  
    Slab_M(const Slab_M& slab);
    ~Slab_M() {}

    Slab_M& operator= (const Slab_M& slab);
    bool  operator==(const Waveguide& w) const;

    std::vector<Material*> get_materials() const
      {return materials;}
    
    bool contains(const Material& m) const
      {return std::find(materials.begin(), materials.end(), &m) 
         != materials.end();}

    bool no_gain_present() const;

    Complex eps_at(const Coord& coord) const;    
    Complex  mu_at(const Coord& coord) const;

    Complex eps_avg() const;

    Complex get_width() const
      {return discontinuities.back();}

    Complex c1_size() const 
      {return discontinuities.back();}

    void find_modes();

    std::vector<Complex> get_params() const;
    void set_params(const std::vector<Complex>&);

  protected:

    std::vector<Complex> find_kt(std::vector<Complex>& old_kt);
    std::vector<Complex> find_kt_from_scratch_by_ADR();
    std::vector<Complex> find_kt_from_scratch_by_track();
    std::vector<Complex> find_kt_by_sweep(std::vector<Complex>& old_kt);
    std::vector<Complex> find_kt_from_estimates();

    cVector estimate_kz2_from_RCWA();
    cVector estimate_kz2_from_uniform_modes();      

    void build_modeset(const std::vector<Complex>& kt);

    std::vector<Complex> params; // Last parameters of dispersion relation.
    
    std::vector<Material*> materials;
    std::vector<Complex>   thicknesses;

    int M_series;

    cVector fourier_eps(int M, cVector* inv_eps=NULL) const;

    friend class Slab_M_Mode;
    friend class UniformSlab;
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
    
    std::vector<Material*> get_materials() const
      {std::vector<Material*> m; m.push_back(core); return m;}

    Complex eps_avg() const 
      {return core->eps();}
    
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

    Complex c1_size() const 
      {return discontinuities[0];}

    void find_modes();

    std::vector<Complex> get_params() const;
    void set_params(const std::vector<Complex>&);

  protected:

    std::vector<Complex> find_kt();
    void build_modeset(std::vector<Complex>& kt);
};



#endif
