
/////////////////////////////////////////////////////////////////////////////
//
// File:     generalslab.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
//
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef GENERALSLAB_H
#define GENERALSLAB_H

#include "../../waveguide.h"
#include "../../expression.h"
#include "../../math/calculus/function.h"

/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: SlabGlobal
//
//   Groups global variables related to slab structures.
//
/////////////////////////////////////////////////////////////////////////////

class SlabWall; // forward declaration - see slabwall.h
struct SlabGlobal
{
    Real      lower_PML;
    Real      upper_PML;
    SlabWall* lowerwall; // NULL: electric wall.
    SlabWall* upperwall;
};

extern SlabGlobal global_slab;



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabImpl
//
//  Class containing the implementation shared by the different types of
//  slabs.
//  
/////////////////////////////////////////////////////////////////////////////

struct SlabCache; // forward declarations - see slaboverlap.h and slabmode.h
class SlabMode;

class SlabImpl : public MultiWaveguide
{ 
  public:

    SlabImpl() : lowerwall(NULL), upperwall(NULL), dummy(0) {}
    ~SlabImpl();

    void set_lower_wall(SlabWall& lower) {lowerwall=&lower;}
    void set_upper_wall(SlabWall& upper) {upperwall=&upper;}
    
    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const;
    
    virtual Complex get_width() const = 0;

    virtual Complex eps_avg() const = 0;
    
    void calc_overlap_matrices
      (MultiWaveguide*, cMatrix*, cMatrix*,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);

    virtual std::vector<Complex> get_params() const = 0;
    virtual void set_params(const std::vector<Complex>&) = 0;

    cVector expand_field(ComplexFunction* f, Real eps);

    std::vector<Complex> disc_intersect(SlabImpl* medium_II);

    void fill_field_cache(SlabCache* cache, SlabImpl* medium_II,
                          const std::vector<Complex>& disc);

    std::vector<Complex> get_discontinuities() const 
      {return discontinuities;}

    void set_dummy(bool b) {dummy = b;}
    bool is_dummy() const {return dummy;}

  protected:

    bool dummy;

    SlabWall* lowerwall; // NULL means use wall from global_slab.
    SlabWall* upperwall;

    // x-values of the interfaces, not including lower wall, including
    // upper wall.
    
    std::vector<Complex> discontinuities;

    friend class Slab_M_Mode;
    friend class UniformSlabMode;

    friend Complex overlap(const SlabMode*, const SlabMode*,
                           const SlabCache* c,
                           const std::vector<Complex>* v,
                           int, int, int, int);

    friend void overlap_TM_TE(const SlabMode*, const SlabMode*,
			      Complex* Ex_Hz, Complex* Ez_Hx,
			      const SlabCache* cache,
			      const std::vector<Complex>* disc,
			      int, int, int, int);
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Slab
//
//   Encapsulates different kinds of slabs through a pointer.
//   Can be initialised by an expression.
//   
/////////////////////////////////////////////////////////////////////////////

class Slab : public MultiWaveguide
{
  public:

    Slab(const Expression& ex);
    Slab(const Term& t);    
    ~Slab() {delete s;}

    bool operator==(const Waveguide& w)    const {return this == &w;}
    std::vector<Material*> get_materials() const {return s->get_materials();}
    bool contains(const Material& m)       const {return s->contains(m);}
    bool no_gain_present()                 const {return s->no_gain_present();}

    void find_modes() {return s->find_modes();}
    
    Mode* get_mode(int i)    const {return s->get_mode(i);}
    Mode* get_fw_mode(int i) const {return s->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return s->get_bw_mode(i);}
    
    void set_lower_wall(SlabWall& lower) const {s->set_lower_wall(lower);}
    void set_upper_wall(SlabWall& upper) const {s->set_upper_wall(upper);}
    
    Complex eps_at(const Coord& coord) const {return s->eps_at(coord);}
    Complex  mu_at(const Coord& coord) const {return s-> mu_at(coord);}

    Complex eps_avg() const {return s->eps_avg();}

    int N() const {return s->N();}

    Complex get_width() const {return s->get_width();}
    Complex   c1_size() const {return s->c1_size();}

    std::vector<Complex> get_params() const {return s->get_params();}
    void set_params(const std::vector<Complex>& p) {s->set_params(p);}

    std::vector<Complex> get_discontinuities() const 
      {return s->get_discontinuities();}

    const FieldExpansion field_from_source
      (const Coord& pos, const Coord& orientation)
        {return s->field_from_source(pos, orientation);}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const
      {return s->S_flux(f, c1_start, c1_stop, precision);}

    void calc_overlap_matrices
      (MultiWaveguide* w2, cMatrix* O_I_II, cMatrix* O_II_I,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL)
      {return s->calc_overlap_matrices(dynamic_cast<Slab*>(w2)->s,
                                       O_I_II,O_II_I,O_I_I,O_II_II);}

    cVector expand_field(ComplexFunction* f, Real eps=1e-4)
      {return s->expand_field(f, eps);}

    SlabImpl* get_impl() const {return s;}
    
    std::string repr() const {return s->repr();}

    void set_dummy(bool b) {s->set_dummy(b);}
    bool is_dummy() const {return s->is_dummy();}
    
  protected:

    SlabImpl* s;
};

inline std::ostream& operator<<(std::ostream& s, const Slab& slab)
  {return s << slab.repr();}



#endif
