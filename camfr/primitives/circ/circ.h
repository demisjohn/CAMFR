
/////////////////////////////////////////////////////////////////////////////
//
// File:     circ.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990513
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CIRC_H
#define CIRC_H

#include "../../util/index.h"
#include "../../waveguide.h"
#include "../../material.h"
#include "../../coord.h"
#include "../../expression.h"
#include "circoverlap.h"
#include "circdisp.h"

/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: CircGlobal
//
//   Groups global variables related to circular structures.
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {cos_type, sin_type} Fieldtype;

struct CircGlobal
{
    Real PML;
    int order; 
    Fieldtype fieldtype;
};

extern CircGlobal global_circ;



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_M
//
//   Cross section consisting of arbitrary number M of rings, surrounded
//   by metal cylinder:
//
//      ring 0   :     0    -->  rad[0]  :  mat[0]
//      ring 1   :  rad[0]  -->  rad[1]  :  mat[1]
//      ..
//      ring M-1 : rad[M-2] --> rad[M-1] : mat[M-1]
//  
/////////////////////////////////////////////////////////////////////////////

// forward declaration - see circ.cpp and circmode.h

struct CircCache;
class Circ_M_Mode;
class Circ_2_Mode;
class Circ_1_Mode;

class Circ_M : public MultiWaveguide
{
  public:

    Circ_M() {}
      
    Circ_M(const std::vector<Complex> &r, const std::vector<Material*> &m)
      : M(m.size()), radius(r), material(m) {}

    Complex c1_size() const {return radius.back();}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const;

    bool operator==(const Waveguide& w) const;

    Complex lateral_S_corr(const Coord& c) const 
      {return (global_circ.order == 0) ? 2*pi*c.c1 : pi*c.c1;} 
    
    bool contains(const Material& m) const;

    bool no_gain_present() const;

    std::vector<Material*> get_materials() const
      {return material;}

    Material* material_at(const Coord& coord) const;

    Complex kt_to_kz(const Complex& kt);
    Complex kz_to_kt(const Complex& kz);
    
    void find_modes();
    
    void calc_overlap_matrices
      (MultiWaveguide* w, cMatrix* O_I_II, cMatrix* O_II_I,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL);
    
    // Create uniform stretching, starting halfway in outer region
    // and giving this particular complex R.
    
    CoordStretcher createStretcher();

  protected:

    void find_modes_from_scratch_by_track();

    static Hankel hankel;

    unsigned int M;
    
    std::vector<Complex>   radius;
    std::vector<Material*> material;

    std::vector<Complex> params;

    friend class CircMode;
    friend class Circ_M_Mode;
    friend class Circ_2_Mode;
    friend class Circ_1_Mode;
    friend Complex overlap(const CircMode*, const CircMode*,
                           const CircCache* c,
                           const std::vector<Complex>* v,
                           int, int, int, int);
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_2
//
//   Core material with radius r,
//   cladding material enclosed in metal cilinder of radius R > r.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_2 : public Circ_M
{
  public:

    Circ_2(const Complex& r, Material& core,
           const Complex& R, Material& cladding);
    
    void find_modes();

  protected:

    void find_modes_from_scratch_by_ADR();
    void find_modes_from_scratch_by_track();
    void find_modes_by_sweep();
    
    std::vector<Complex> kr2_backward;
    std::vector<Complex> guided_disp_params;
    std::vector<Complex> rad_disp_params;

    unsigned int no_of_guided_modes;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_1
//
//   Uniform material enclosed in metal cilinder of radius R.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_1 : public Circ_M
{ 
  public:

    Circ_1(const Complex& radius, Material& m);

    void find_modes();
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ
//
//   Encapsulates different kinds of Circ structures through a pointer.
//   Can be initialised by an expression.
//   
/////////////////////////////////////////////////////////////////////////////

class Circ : public MultiWaveguide
{
  public:

    Circ(const Expression& ex);
    Circ(const Term& t);    
    ~Circ() {delete c;}

    bool operator==(const Waveguide& w)    const {return this == &w;}
    std::vector<Material*> get_materials() const {return c->get_materials();}
    bool contains(const Material& m)       const {return c->contains(m);}
    bool no_gain_present()                 const {return c->no_gain_present();}
    
    Complex c1_size() const {return c->c1_size();}

    void find_modes() {return c->find_modes();}
    
    Mode* get_mode(int i)    const {return c->get_mode(i);}
    Mode* get_fw_mode(int i) const {return c->get_fw_mode(i);}
    Mode* get_bw_mode(int i) const {return c->get_bw_mode(i);}
    
    Material* material_at(const Coord& coord) const 
      {return c->material_at(coord);}
    
    int N() const {return c->N();}

    const FieldExpansion field_from_source
      (const Coord& pos, const Coord& orientation)
        {return c->field_from_source(pos, orientation);}

    Real S_flux(const FieldExpansion& f,
                Real c1_start, Real c1_stop,
                Real precision = 1e-10) const
      {return c->S_flux(f, c1_start, c1_stop, precision);}

    void calc_overlap_matrices
      (MultiWaveguide* w2, cMatrix* O_I_II, cMatrix* O_II_I,
       cMatrix* O_I_I=NULL, cMatrix* O_II_II=NULL)
      {return c->calc_overlap_matrices(dynamic_cast<Circ*>(w2)->c,
                                       O_I_II,O_II_I,O_I_I,O_II_II);}

    std::string repr() const {return c->repr();}
    
  protected:

    Circ_M* c;
};

inline std::ostream& operator<<(std::ostream& s, const Circ& circ)
  {return s << circ.repr();}  



#endif
