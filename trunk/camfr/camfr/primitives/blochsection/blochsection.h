
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsection.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20050517
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef BLOCHSECTION_H
#define BLOCHSECTION_H

#include <string>
#include <vector>
#include "../../stack.h"
#include "../slab/generalslab.h"


/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: BlochSectionGlobal
//
//   Groups global variables related to section structures.
//
/////////////////////////////////////////////////////////////////////////////

// TODO: break out

struct BlochSectionGlobal
{
    int  Mx;
    int  My;
    Complex alpha0;
    Complex beta0;
    Real left_PML;
    Real right_PML;
    Real PML_fraction;
};

extern BlochSectionGlobal global_blochsection;



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSectionImpl
//
//  Class containing the implementation shared by the different types of
//  bloch sections.
//  
/////////////////////////////////////////////////////////////////////////////

class BlochSectionMode;

class BlochSectionImpl : public MultiWaveguide
{ 
  public:

    BlochSectionImpl() {}

    virtual Complex get_width()  const = 0;
    virtual Complex get_height() const = 0;
    virtual Complex c1_size()    const {return get_width();}

    std::vector<Complex> get_disc() const {return discontinuities;}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const
      {std::cout << "Not yet implemented." << std::endl; return -99999;}

    void calc_overlap_matrices
      (MultiWaveguide*, cMatrix*, cMatrix*,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);

  protected:

    // z-values of the interfaces, not including left wall, including
    // right wall.
    
    std::vector<Complex> discontinuities;
    std::vector<Slab*> slabs;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSection
//
//   A 2D periodic grating with an arbitrary 2D cartesian cross section.
//   Encapsulates different kinds of Blochsections through a pointer.
//   
/////////////////////////////////////////////////////////////////////////////

class BlochSection : public MultiWaveguide
{
  public:

    BlochSection(const Term& t);
    BlochSection(Expression& ex);
    ~BlochSection() {delete s;}

    bool operator==(const Waveguide& w)    const {return this == &w;}
    std::vector<Material*> get_materials() const {return s->get_materials();}
    bool contains(const Material& m)       const {return s->contains(m);}
    bool no_gain_present()                 const {return s->no_gain_present();}

    void find_modes() {return s->find_modes();}
    
    Mode* get_mode(int i)    const {return s->get_mode(i);}
    Mode* get_fw_mode(int i) const {return s->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return s->get_bw_mode(i);}
    
    Material* material_at(const Coord& co) const {return s->material_at(co);}
    
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
      {return s->calc_overlap_matrices(dynamic_cast<BlochSection*>(w2)->s,
                                       O_I_II,O_II_I,O_I_I,O_II_II);}
    
    std::string repr() const {return s->repr();}
    
  protected:

    BlochSectionImpl* s;
};

inline std::ostream& operator<<(std::ostream& s, const BlochSection& section)
  {return s << section.repr();}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSection2D
//
//   Most general section: 2D cartesion section consisting of a horizontal 
//   sequence in the x-direction of several Slabs. Each slab is oriented 
//   vertically along y.
//  
/////////////////////////////////////////////////////////////////////////////

class BlochSection2D : public BlochSectionImpl
{ 
  public:

    BlochSection2D() {}
    BlochSection2D(Expression& ex);   
    BlochSection2D(const BlochSection2D& section);

    BlochSection2D& operator= (const BlochSection2D& section);

    bool operator== (const Waveguide& w) const
      {return (this == &w);}

    std::vector<Material*> get_materials() const 
      {return materials;}    

    bool contains(const Material& m) const
      {return std::find(materials.begin(), materials.end(), &m) 
         != materials.end();}

    bool no_gain_present() const;

    Material* material_at(const Coord& coord) const;

    Complex get_width() const
      {return st.get_total_thickness();}

    Complex get_height() const
      {return slabs[0]->get_width();}

    void find_modes();

  protected:

    // TODO: see if we can remove st.

    Stack st;

    void create_FG_li(cMatrix* F, cMatrix* G, int M, int N,
                      const Complex& alpha0, const Complex& beta0);    

    void create_FG_li_biaxial(cMatrix* F, cMatrix* G, int M, int N,
                      const Complex& alpha0, const Complex& beta0);

    std::vector<Material*> materials;

    friend class BlochSection2D_Mode;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: UniformBlochSection
//  
/////////////////////////////////////////////////////////////////////////////

class UniformBlochSection : public BlochSectionImpl
{ 
  public:

    UniformBlochSection(Material* mat, const Complex& W_, const Complex& H_) 
      : W(W_), H(H_) {core = mat;}
    
    bool operator== (const Waveguide& w) const
      {return this == &w;}
    
    std::vector<Material*> get_materials() const
      {std::vector<Material*> m; m.push_back(core); return m;}
    
    bool contains(const Material& m) const
      {return m == *core;}

    bool no_gain_present() const
      {return core->no_gain_present();}

    Material* material_at(const Coord& c) const
      {return core;}

    Complex get_width() const
      {return W;}

    Complex get_height() const
      {return H;}

    void find_modes();

  protected:

    Complex W, H;
};



#endif
