
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
#include "../../stack.h"
#include "../slab/generalslab.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SectionImpl
//
//  Class containing the implementation shared by the different types of
//  sections.
//  
/////////////////////////////////////////////////////////////////////////////

class SectionImpl : public MultiWaveguide
{ 
  public:

    SectionImpl() {}
     
    virtual Complex get_width()  const = 0;
    virtual Complex get_height() const = 0;
    virtual Complex c1_size()    const {return get_width();}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const
      {std::cout << "Not yet implemented." << std::endl; return -99999;}

    void calc_overlap_matrices
      (MultiWaveguide*, cMatrix*, cMatrix*,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);
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

    Section(const Expression& ex, int M=global.N);
    Section(const Expression& left_ex, const Expression& right_ex, 
            int M=global.N);
    Section(const Term& t);    
    ~Section() {delete s;}

    bool operator==(const Waveguide& w)    const {return this == &w;}
    std::vector<Material*> get_materials() const {return s->get_materials();}
    bool contains(const Material& m)       const {return s->contains(m);}
    bool no_gain_present()                 const {return s->no_gain_present();}

    void find_modes() {return s->find_modes();}
    
    Mode* get_mode(int i)    const {return s->get_mode(i);}
    Mode* get_fw_mode(int i) const {return s->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return s->get_bw_mode(i);}
    
    Complex eps_at(const Coord& coord) const {return s->eps_at(coord);}
    Complex  mu_at(const Coord& coord) const {return s-> mu_at(coord);}
    
    int N() const {return s->N();}

    Complex get_width()  const {return s->get_width();}
    Complex get_height() const {return s->get_height();}
    Complex c1_size()    const {return s->c1_size();}

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
    Section2D(const Expression& left_ex, const Expression& right_ex, int M);   
    Section2D(const Section2D& section);
    ~Section2D() {}

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
      {return dynamic_cast<Slab*>(left.get_inc())->get_width();}

    int get_M() const {return M;}
    
    void find_modes();

  protected:

    int M;

    Stack left;
    Stack right;

    bool symmetric;

    void find_modes_from_scratch_by_ADR();
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

    Section1D(Slab& slab, const Complex& width) 
      : s(&slab), d(width) {uniform = false; core = s->get_core();}

    ~Section1D() {}
    
    bool operator== (const Waveguide& w) const
      {return this == &w;}
    
    std::vector<Material*> get_materials() const
      {return s->get_materials();}
    
    bool contains(const Material& m) const
      {return s->contains(m);}

    bool no_gain_present() const
      {return s->no_gain_present();}

    Complex eps_at(const Coord& c) const
      {return s->eps_at(Coord(c.c2, 0.0, 0.0, c.c2_limit, Plus, Plus));}

    Complex  mu_at(const Coord& c) const
      {return s->mu_at (Coord(c.c2, 0.0, 0.0, c.c2_limit, Plus, Plus));}

    Complex get_width() const
      {return d;}

    Complex get_height() const
      {return s->get_width();}

    void find_modes();

  protected:

    Slab* s;
    Complex d;
};



#endif
