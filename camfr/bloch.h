
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

class BlochStack : public MultiWaveguide
{
  public:

    BlochStack(const Expression& e);
    
    Complex eps_at(const Coord& coord) const {return stack.eps_at(coord);}
    Complex  mu_at(const Coord& coord) const {return stack. mu_at(coord);}

    bool operator==(const Waveguide& w) const {return &w==this;}
    vector<Material*>   get_materials() const {return stack.get_materials();}
    bool    contains(const Material& m) const {return stack.contains(m);}
    bool              no_gain_present() const {return stack.no_gain_present();}

    Complex get_total_thickness() const {return stack.get_total_thickness();}

    void find_modes();

    cVector get_beta_vector() const;
    
  protected:

    Stack stack;

  private:

    void find_modes_GEV();
    void find_modes_T();

    void find_modes_diag();

    void brillouin_eliminate_modes();
};



/////////////////////////////////////////////////////////////////////////////
//
// BlochMode
//
/////////////////////////////////////////////////////////////////////////////

class BlochMode : public Mode
{
  public:

    BlochMode(const Polarisation pol, const Complex& kz, Stack* s,
              cVector& F, cVector& B);

    Field field(const Coord& coord) const;

    cVector fw_field() const {return interface_field[0].fw;}
    cVector bw_field() const {return interface_field[0].bw;}    

    Real S_flux(Real c1_start, Real c1_stop, Real eps=1e-10) const
      {return dynamic_cast<MultiWaveguide*>(geom->get_inc())->
         S_flux(interface_field[0], c1_start, c1_stop, eps);}

  protected:

    Stack* geom;

    mutable vector<FieldExpansion> interface_field;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS Bandtracer
//
/////////////////////////////////////////////////////////////////////////////

class Bandtracer : public Tracesorter
{
  public:

    Bandtracer::Bandtracer(Real period) 
      {
        add_turning_point( pi/period);
        add_turning_point(-pi/period);
      }
};



#endif



