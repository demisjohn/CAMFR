
/////////////////////////////////////////////////////////////////////////////
//
// File:     section.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020129
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SECTION_H
#define SECTION_H

#include <string>
#include <vector>
#include "../../stack.h"
#include "../slab/generalslab.h"

#include "sectiondisp.h" // TMP

/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: SectionGlobal
//
//   Groups global variables related to section structures.
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {E_wall, H_wall} Section_wall_type;

struct SectionGlobal
{
    Real  left_PML;
    Real right_PML;
    Section_wall_type  leftwall;
    Section_wall_type rightwall;
    bool guided_only; // TMP variable?
};

extern SectionGlobal global_section;



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SectionImpl
//
//  Class containing the implementation shared by the different types of
//  sections.
//  
/////////////////////////////////////////////////////////////////////////////

class SectionMode;
class SectionCache;
class OverlapMatrices;
class RefSection;

class SectionImpl : public MultiWaveguide
{ 
  public:

    SectionImpl() {}

    virtual Complex get_width()  const = 0;
    virtual Complex get_height() const = 0;
    virtual Complex c1_size()    const {return get_width();}

    int get_M1() const {return M1;}
    int get_M2() const {return M2;}

    std::vector<Complex> get_disc() const {return discontinuities;}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const
      {std::cout << "Not yet implemented." << std::endl; return -99999;}

    void calc_overlap_matrices
      (MultiWaveguide*, cMatrix*, cMatrix*,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);

    SectionDisp* get_disp() const {return f;} // TMP

  protected:
    
  public: // tmp

    int M1;
    int M2;

    SectionDisp* f;

    // z-values of the interfaces, not including left wall, including
    // right wall.
    
    std::vector<Complex> discontinuities;
    std::vector<Slab*> slabs;

    friend Complex overlap_slice(SectionMode*, SectionMode*,
                                 const Complex&, const Complex&,
                                 FieldExpansion*, FieldExpansion*, 
                                 OverlapMatrices*, int, int);

    friend class RefSection;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section
//
//   A waveguide with an arbitrary 2D cartesian cross section.
//   Encapsulates different kinds of sections through a pointer.
//   
/////////////////////////////////////////////////////////////////////////////

class Section : public MultiWaveguide
{
  public:

    Section(const Term& t);
    Section(Expression& ex, 
            int M1=int(global.N*global.mode_surplus), int M2=global.N);
    Section(Expression& left_ex, Expression& right_ex, 
            int M1=int(global.N*global.mode_surplus), int M2=global.N);
    ~Section() {delete s; delete leftwall_sc; delete rightwall_sc;}

    bool operator==(const Waveguide& w)    const {return this == &w;}
    std::vector<Material*> get_materials() const {return s->get_materials();}
    bool contains(const Material& m)       const {return s->contains(m);}
    bool no_gain_present()                 const {return s->no_gain_present();}

    void find_modes() {return s->find_modes();}
    
    Mode* get_mode(int i)    const {return s->get_mode(i);}
    Mode* get_fw_mode(int i) const {return s->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return s->get_bw_mode(i);}
    
    Complex eps_at(const Coord& co) const {return s->eps_at(co);}
    Complex  mu_at(const Coord& co) const {return s-> mu_at(co);}
    Complex   n_at(const Coord& co) const {return sqrt(eps_at(co)/eps0);}  
    
    int N() const {return s->N();}

    Complex get_width()  const {return s->get_width();}
    Complex get_height() const {return s->get_height();}
    Complex c1_size()    const {return s->c1_size();}

    Complex get_disp(const Complex& z) const 
      {return (*s->get_disp())(z);} // TMP

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
      {return s->calc_overlap_matrices(dynamic_cast<Section*>(w2)->s,
                                       O_I_II,O_II_I,O_I_I,O_II_II);}
    
    std::string repr() const {return s->repr();}
    
  protected:

    SectionImpl* s;
    DiagScatterer* leftwall_sc;
    DiagScatterer* rightwall_sc;
};

inline std::ostream& operator<<(std::ostream& s, const Section& section)
  {return s << section.repr();}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section2D
//
//   Most general section: 2D cartesion section consisting of a horizontal 
//   sequence in the x-direction of several Slabs. Each slab is oriented 
//   vertically along y.
//  
/////////////////////////////////////////////////////////////////////////////

class Section2D : public SectionImpl
{ 
  public:

    Section2D() {}
    Section2D(Expression& left_ex, Expression& right_ex, int M1, int M2);   
    Section2D(const Section2D& section);

    Section2D& operator= (const Section2D& section);

    bool operator== (const Waveguide& w) const
      {return (this == &w);}

    std::vector<Material*> get_materials() const 
      {return materials;}    

    bool contains(const Material& m) const
      {return std::find(materials.begin(), materials.end(), &m) 
         != materials.end();}

    bool no_gain_present() const;

    Complex eps_at(const Coord& coord) const;
    Complex  mu_at(const Coord& coord) const;

    Complex get_width() const
      {return left.get_total_thickness() + right.get_total_thickness();}

    Complex get_height() const
      {return slabs[0]->get_width();}
    
    void find_modes();

  protected:

    Stack left;
    Stack right;

    bool symmetric;

    void find_modes_from_series();
    void find_modes_from_scratch_by_track();
    void find_modes_by_sweep();

    std::vector<Complex> params; // Last parameters of dispersion relation.
    std::vector<Material*> materials;

    friend class Section2D_Mode;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section1D
//
//   Special case: 2D cartesian section where the index profile only varies 
//   in the vertical y-direction (i.e. stack of uniform horizontal layers).
//  
/////////////////////////////////////////////////////////////////////////////

class Section1D : public SectionImpl
{ 
  public:

    Section1D(Slab& slab_, const Complex& width) 
      : slab(&slab_), d(width) 
        {uniform=false; core=slab->get_core(); M1=0; M2=1;}

    ~Section1D() {}
    
    bool operator== (const Waveguide& w) const
      {return this == &w;}
    
    std::vector<Material*> get_materials() const
      {return slab->get_materials();}
    
    bool contains(const Material& m) const
      {return slab->contains(m);}

    bool no_gain_present() const
      {return slab->no_gain_present();}

    Complex eps_at(const Coord& c) const
      {return slab->eps_at(Coord(c.c2, 0.0, 0.0, c.c2_limit, Plus, Plus));}

    Complex  mu_at(const Coord& c) const
      {return slab->mu_at (Coord(c.c2, 0.0, 0.0, c.c2_limit, Plus, Plus));}

    Complex get_width() const
      {return d;}

    Complex get_height() const
      {return slab->get_width();}

    void find_modes();

  protected:

    Slab* slab;
    Complex d;
};



#endif
