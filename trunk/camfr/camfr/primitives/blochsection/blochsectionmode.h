
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsectionmode.h
// Author:   Peter.Bienstman@UGent.be
// Date:     20050518
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef BLOCHSECTIONMODE_H
#define BLOCHSECTIONMODE_H

#include <vector>
#include "../../mode.h"
#include "blochsection.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSectionMode
//
//   Interface class for BlochSection2D_Mode and BlochSection1D_Mode.
//
/////////////////////////////////////////////////////////////////////////////

class BlochSectionMode : public Mode
{
  public:

    BlochSectionMode(Polarisation pol, const Complex& kz0, 
                     BlochSectionImpl* geom_)
      : Mode(pol, kz0, -kz0), geom(geom_) {}

    BlochSectionImpl* get_geom() const {return geom;}

    virtual void normalise() = 0;

  protected:

    BlochSectionImpl* geom;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSection2D_Mode
//
/////////////////////////////////////////////////////////////////////////////

class BlochSection2D_Mode : public BlochSectionMode
{
  public:

    BlochSection2D_Mode(Polarisation pol, const Complex& kz, 
                        BlochSectionImpl* geom,
                        const cVector& Ex, const cVector& Ey,
                        const cVector& Hx, const cVector& Hy);

    Field field(const Coord& coord) const;
    
    void normalise();

    friend Complex overlap(const BlochSection2D_Mode* sec_I_mode,
                           const BlochSection2D_Mode* sec_II_mode);

  protected:

    cVector Ex, Ey, Hx, Hy;
};


// TODO: Implement.

#if 0

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: BlochSection1D_Mode
//
/////////////////////////////////////////////////////////////////////////////

class BlochSection1D_Mode : public BlochSectionMode
{
  public:

    BlochSection1D_Mode(Polarisation pol, const Complex& kz, 
                   SlabMode* m, BlochSection1D* geom);

    Field field(const Coord& coord) const;

    void normalise();

  protected:

    SlabMode* m;
    Complex fw0, bw0;
};
#endif



#endif



