
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabmode.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000203
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SLABMODE_H
#define SLABMODE_H

#include <vector>
#include "../../../mode.h"
#include "slab.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabMode
//
//   Interface class for Slab_M_Mode and UniformSlabMode
//
/////////////////////////////////////////////////////////////////////////////

class SlabMode : public Mode
{
  public:

    SlabMode(Polarisation pol, const Complex& kz0, const SlabImpl* geom_)
      : Mode(pol, kz0, -kz0), geom(geom_) {}

    Complex get_kz0() const {return kz;}
    Complex get_kz () const;

    Complex get_sin() const {return global.slab_ky / kz;}
    Complex get_cos() const;

    virtual void normalise() = 0;

    const SlabImpl* get_geom() const {return geom;}

    // Calculate amplitudes of forward(right) and backward(left) waves.
    
    virtual void forw_backw_at(const Coord& coord,
                               Complex* fw, Complex* bw) const = 0;

    // Calculate total (forward+backward) field.
    
    Field field(const Coord& coord) const;

    virtual Complex kx_at(const Coord& coord) const = 0;

  protected:

    const SlabImpl* geom;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Slab_M_Mode
//
/////////////////////////////////////////////////////////////////////////////

class Slab_M_Mode : public SlabMode
{
  public:

    Slab_M_Mode(Polarisation pol,  const Complex& kz, 
                const Complex& kt, const Slab_M* geom,
                bool calc_fw = true);

    Slab_M_Mode(Polarisation pol,   const Complex& kz,
                const Complex& kt,  const Slab_M* geom,
                Stack& lower_stack, Stack& upper_stack);   

    void normalise();
    
    void forw_backw_at(const Coord& coord, Complex* fw, Complex* bw) const;

    Complex kx_at(const Coord& c) const
      {return kx[index_lookup(c.c1, c.c1_limit, geom->discontinuities)];}

    Complex get_kt() const {return kt;}
    
  protected:

    void calc_fw_bw(bool calc_fw = true);
    void calc_fw_bw(Stack& lower_stack, Stack& upper_stack);

    // Expansion coeff. of forward and backward waves at each interface.
        
    std::vector<Complex> fw; // towards x=+inf
    std::vector<Complex> bw; // towards x=-inf

    // x-component of wavevector.
    
    std::vector<Complex> kx;

    // Root of the orignal dispersion relation.
    
    Complex kt;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: UniformSlabMode
//
/////////////////////////////////////////////////////////////////////////////

class UniformSlabMode : public SlabMode
{
  public:

    UniformSlabMode(Polarisation pol, const Complex& kz,
                    const UniformSlab* geom);

    void normalise();
    
    void forw_backw_at(const Coord& coord, Complex* fw, Complex* bw) const;

    Complex kx_at(const Coord& coord) const
      {return kx;}

  protected:

    Complex kx;
    Complex fw0;
    Complex bw0;
};



#endif



