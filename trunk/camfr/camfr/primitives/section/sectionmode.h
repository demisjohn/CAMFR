
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

    SectionMode(Polarisation pol, const Complex& kz0, SectionImpl* geom_)
      : Mode(pol, kz0, -kz0), geom(geom_) {}

    SectionImpl* get_geom() const {return geom;}

    virtual void get_fw_bw(const Complex& c, Limit c_limit,
                           cVector* fw, cVector* bw) const = 0;

    virtual void normalise() = 0;

  protected:

    SectionImpl* geom;

};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section2D_Mode
//
/////////////////////////////////////////////////////////////////////////////

class Section2D_Mode : public SectionMode
{
  public:

    Section2D_Mode(Polarisation pol, const Complex& kz, Section2D* geom,
                   cVector* Ex=0,cVector* Ey=0,cVector* Hx=0,cVector* Hy=0,
                   bool corrected=true);
  
    ~Section2D_Mode();

    Field field(const Coord& coord) const;

    void get_fw_bw(const Complex&, Limit, cVector*, cVector*) const;
    
    void normalise();

    friend Complex overlap_pw(const Section2D_Mode* sec_I_mode, 
                              const Section2D_Mode* sec_II_mode);    

    friend Complex overlap(const Section2D_Mode* sec_I_mode,
                           const Section2D_Mode* sec_II_mode);

  protected:

    mutable std::vector<FieldExpansion>  left_interface_field;
    mutable std::vector<FieldExpansion> right_interface_field;

    cVector *Ex, *Ey, *Hx, *Hy;

    bool corrected;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Section1D_Mode
//
/////////////////////////////////////////////////////////////////////////////

class Section1D_Mode : public SectionMode
{
  public:

    Section1D_Mode(Polarisation pol, const Complex& kz, 
                   SlabMode* m, Section1D* geom);

    Field field(const Coord& coord) const;

    void get_fw_bw(const Complex&, Limit, cVector*, cVector*) const;

    void normalise();

  protected:

    void get_fw_bw(const Complex& c, Complex* fw, Complex* bw) const;

    SlabMode* m;
    Complex fw0, bw0;
};



#endif



