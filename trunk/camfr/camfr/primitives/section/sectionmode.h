
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectionmode.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020225
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SECTIONMODE_H
#define SECTIONMODE_H

#include <vector>
#include "../../mode.h"
#include "section.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SectionMode
//
//   Interface class for Section2D_Mode and Section1D_Mode.
//
/////////////////////////////////////////////////////////////////////////////

class SectionMode : public Mode
{
  public:

    SectionMode(Polarisation pol,const Complex& kz0,const SectionImpl* geom_)
      : Mode(pol, kz0, -kz0), geom(geom_) {}

    Complex get_kz () const {return kz;}

    const SectionImpl* get_geom() const {return geom;}
    
    virtual Field field(const Coord& coord) const;

    virtual void normalise() = 0;

  protected:

    const SectionImpl* geom;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section2D_Mode
//
/////////////////////////////////////////////////////////////////////////////

class Section2D_Mode : public SectionMode
{
  public:

    Section2D_Mode(Polarisation pol, const Complex& kz, 
                   const Section2D* geom, const cVector& fw_field);

    Field field(const Coord& coord) const;
    
    void normalise();

  protected:

    mutable vector<FieldExpansion>  left_interface_field;
    mutable vector<FieldExpansion> right_interface_field;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: UniformSectionMode
//
/////////////////////////////////////////////////////////////////////////////

class Section1D_Mode : public SectionMode
{
  public:

    Section1D_Mode(Polarisation pol, const Complex& kz,
                   const Section1D* geom);

    Field field(const Coord& coord) const;

    void normalise();
};



#endif



