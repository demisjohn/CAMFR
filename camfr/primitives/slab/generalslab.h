
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
    SlabWall*  leftwall; // NULL: electric wall.
    SlabWall* rightwall;
    Complex beta; // The out-of plane component for off-angle
                  // propagation k_y, not k_z.
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

    SlabImpl() : leftwall(NULL), rightwall(NULL) {}
    ~SlabImpl();

    void  set_left_wall(SlabWall&  left) { leftwall=&left;}
    void set_right_wall(SlabWall& right) {rightwall=&right;}
    
    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const;
    
    virtual Complex get_width() const = 0;
    
    void calc_overlap_matrices
      (MultiWaveguide*, cMatrix*, cMatrix*,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);

    virtual vector<Complex> get_params() const = 0;
    virtual void set_params(const vector<Complex>&) = 0;

  protected:

    SlabWall*  leftwall; // NULL means use wall from global_slab.
    SlabWall* rightwall;

    // x-values of the interfaces, not including left wall, including
    // right wall.
    
    vector<Complex> discontinuities;

    friend class Slab_M_Mode;
    friend class UniformSlabMode;

    friend Complex overlap(const SlabMode*, const SlabMode*,
                           const SlabCache* c,
                           const vector<Complex>* v,
                           int, int, int, int);

    friend void overlap_TM_TE(const SlabMode*, const SlabMode*,
			      Complex* Ex_Hz, Complex* Ez_Hx,
			      const SlabCache* cache,
			      const vector<Complex>* disc,
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

    bool operator==(const Waveguide& w) const {return this == &w;}
    vector<Material*> get_materials()   const {return s->get_materials();}
    bool contains(const Material& m)    const {return s->contains(m);}
    bool no_gain_present()              const {return s->no_gain_present();}

    void find_modes() {return s->find_modes();}
    
    Mode* get_mode(int i)    const {return s->get_mode(i);}
    Mode* get_fw_mode(int i) const {return s->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return s->get_bw_mode(i);}
    
    void  set_left_wall(SlabWall&  left) const {s->set_left_wall(left);}
    void set_right_wall(SlabWall& right) const {s->set_right_wall(right);}
    
    Complex eps_at(const Coord& coord) const {return s->eps_at(coord);}
    Complex  mu_at(const Coord& coord) const {return s-> mu_at(coord);}
    
    int N() const {return s->N();}

    Complex get_width() const {return s->get_width();}

    vector<Complex> get_params() const {return s->get_params();}
    void set_params(const vector<Complex>& p) {s->set_params(p);}

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
    
    string repr() const {return s->repr();}
    
  protected:

    SlabImpl* s;
};

inline ostream& operator<<(ostream& s, const Slab& slab)
  {return s << slab.repr();}



#endif
