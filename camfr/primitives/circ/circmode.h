
/////////////////////////////////////////////////////////////////////////////
//
// File:     circmode.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991104
// Version:  1.2
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CIRCMODE_H
#define CIRCMODE_H

#include "../../mode.h"
#include "circ.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: CircMode
//
//   Mode in a circular symmetric medium surrounded by metal cylinder.
//  
/////////////////////////////////////////////////////////////////////////////

class CircMode : public Mode
{
  public:

    CircMode(Polarisation pol, const Complex& kz, const Circ_M* geom);

    virtual void normalise();

    Field field(const Coord& coord) const {return field_at(coord);}
    
    const Circ_M* get_geom() const {return geom;}
    
    // Useful routines for overlap integrals.
    
    Complex kr_at(const Coord& coord) const
      {return kr[index_lookup(coord.c1, coord.c1_limit, geom->radius)];}
    
    Field field_at(const Coord& coord, Complex* dEzdr=0,
                   Complex* dHzdr=0, bool ang_dep=true) const
      {return field_ring(index_lookup(coord.c1, coord.c1_limit, geom->radius),
                         coord, dEzdr, dHzdr, ang_dep);}

    // ang_dep = false is used by the calculation of the overlap integrals
    // to factor out the phi dependence of the z-components of the fields.
    // Other field components are then incorrect, but are not required by
    // the overlap integral calculation routine.    
    
  protected:
    
    virtual Field field_ring(int i,
                             const Coord& coord,
                             Complex* dEzdr=0,
                             Complex* dHzdr=0,
                             bool ang_dep=true) const = 0;

    const Circ_M* geom;
        
    std::vector<Complex> kr; // Radial component of wavevector in each ring.
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_M_Mode
//
//   Mode in a circular symmetric medium with M rings surrounded by 
//   metal cylinder.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_M_Mode : public CircMode
{

  public:

    Circ_M_Mode(Polarisation pol, const Complex& kz, const Circ_M* geom);

    void normalise();

  protected:

    Field field_ring(int i, const Coord& coord,
                     Complex* dEzdr=0, Complex* dHzdr=0,
                     bool ang_dep=true) const;

    std::vector<cMatrix> T;
    std::vector<cVector> amplitudes;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_2_Mode
//
//   Mode in a circular symmetric medium with core and cladding,
//   surrounded by metal cylinder.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_2_Mode : public CircMode
{
  public:

    Circ_2_Mode(Polarisation pol, const Complex& kz, const Complex& kr1,
                const Complex& kr2, const Circ_M* geom);
    
  protected:
    
    Field field_ring
      (int i, const Coord& coord, Complex* dEzdr=0, Complex* dHzdr=0,
       bool ang_dep=true) const
         {return (i==0) ? field_core    (coord, dEzdr, dHzdr, ang_dep)
                        : field_cladding(coord, dEzdr, dHzdr, ang_dep);}

    Field field_core(const Coord& coord, Complex* dEzdr=0,
                     Complex* dHzdr=0, bool ang_dep=true) const;
    
    Field field_cladding(const Coord& coord, Complex* dEzdr=0,
                         Complex* dHzdr=0,  bool ang_dep=true) const;

    // Exponential scaling of bessel functions in intermediate results.
    
    bool scaling_co; // core
    bool scaling_cl; // cladding 
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Circ_1_Mode
//
//   Mode in a uniform circular symmetric medium
//   surrounded by metal cylinder.
//  
/////////////////////////////////////////////////////////////////////////////

class Circ_1_Mode : public CircMode
{
  public:

    Circ_1_Mode(Polarisation pol,  const Complex& kz,
                const Complex& kr, const Circ_M*  geom);

    void normalise();

  protected:

    Field field_ring(int i, const Coord& coord,
                     Complex* dEzdr=0, Complex* dHzdr=0,
                     bool ang_dep=true) const;
};



#endif



