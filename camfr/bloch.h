
/////////////////////////////////////////////////////////////////////////////
//
// File:     bloch.h
// Authors:  Peter.Bienstman@rug.ac.be, Lieven.Vanholme@rug.ac.be
// Date:     20020207
// Version:  2.0
//
// Copyright (C) 2000-2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef BLOCH_H
#define BLOCH_H

#include "stack.h"
#include "util/tracesorter.h"

/////////////////////////////////////////////////////////////////////////////
//
// BlochStack
//
//   Stack that is repeated an infinite number of times.
//
/////////////////////////////////////////////////////////////////////////////

class BlochMode;

class BlochStack : public MultiWaveguide
{
  public:

    BlochStack(const Expression& e);
    
    Material* material_at(const Coord& coord) const 
      {return stack.material_at(coord);}

    bool operator==(const Waveguide& w) const 
      {return &w==this;}

    std::vector<Material*> get_materials() const
      {return stack.get_materials();}

    bool contains(const Material& m) const 
      {return stack.contains(m);}

    bool no_gain_present() const 
      {return stack.no_gain_present();}

    Complex get_total_thickness() const 
      {return stack.get_total_thickness();}

    Complex c1_size() const 
      {return stack.get_inc()->c1_size();}

    void find_modes();

    cVector get_beta_vector() const;
    
    virtual Mode* get_fw_mode(int i) const;
    virtual Mode* get_bw_mode(int i) const;
    
    // Fool interface building code.

    const Waveguide* get_inc() const {return this;}
    const Waveguide* get_ext() const {return this;}

    // Get true interface waveguide.

    Waveguide* get_inc_wg() const {return stack.get_inc();}
    Waveguide* get_ext_wg() const {return stack.get_ext();} 

    void get_expansion_matrices(cMatrix& ff, cMatrix& fb, 
                                cMatrix& bf, cMatrix& bb, bool left); 
        
  protected:

    Stack stack;

  private:

    void find_modes_GEV();
    void find_modes_T();

    void find_modes_diag();
};



/////////////////////////////////////////////////////////////////////////////
//
// BlochMode
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {undefined,forward,backward} Direction;

class BlochMode : public Mode
{
  public:

    BlochMode(const Polarisation pol, const Complex& kz, Stack* s,
              cVector& F, cVector& B);

    Field field(const Coord& coord) const;

    void fw_bw_field(const Coord& coord, cVector* fw, cVector* bw);

    cVector fw_field() const {return interface_field[0].fw;}
    cVector bw_field() const {return interface_field[0].bw;}    

    Real S_flux(Real c1_start, Real c1_stop, Real eps=1e-10) const
      {return dynamic_cast<MultiWaveguide*>(geom->get_inc())->
         S_flux(interface_field[0], c1_start, c1_stop, eps);}

    Stack* get_geom() const {return geom;}

    Direction get_direction() const {return direction;}
    
  protected:

    Stack* geom;

    mutable std::vector<FieldExpansion> interface_field;

    Direction direction;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS Bandtracer
//
/////////////////////////////////////////////////////////////////////////////

class Bandtracer : public Tracesorter
{
  public:

    Bandtracer(Real period) 
      {
        add_turning_point( pi/period);
        add_turning_point(-pi/period);
      }
};



#endif



