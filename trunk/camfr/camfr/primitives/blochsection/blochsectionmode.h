
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
/////////////////////////////////////////////////////////////////////////////

class BlochSectionMode : public Mode
{
  public:

    BlochSectionMode(Polarisation pol, const Complex& kz, 
                     BlochSectionImpl* geom_,
                     const cVector& Ex, const cVector& Ey,
                     const cVector& Hx, const cVector& Hy);

    BlochSectionImpl* get_geom() const {return geom;}

    Field field(const Coord& coord) const;

    void normalise();

    virtual int get_Mx() const {return -999;}
    virtual int get_My() const {return -999;}

    virtual Complex get_kx() const {return -999;}
    virtual Complex get_ky() const {return -999;}    

    friend Complex overlap(const BlochSectionMode* sec_I_mode,
                           const BlochSectionMode* sec_II_mode);

  protected:

    BlochSectionImpl* geom;    

    // Note: for uniform structures this could be simplified to a scalar,
    // but the gains are probably not very large.

    cVector Ex, Ey, Hx, Hy;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: UniformBlochSectionMode
//
/////////////////////////////////////////////////////////////////////////////

class UniformBlochSectionMode : public BlochSectionMode
{
  public:

    UniformBlochSectionMode(Polarisation pol, const Complex& kz,
                            BlochSectionImpl* geom, int M_, int N_,
                            const cVector& Ex, const cVector& Ey,
                            const cVector& Hx, const cVector& Hy);

    int get_Mx() const {return M;}
    int get_My() const {return N;}    

    Complex get_kx() const;
    Complex get_ky() const;


  protected:
    int M, N;
};



#endif



