
/////////////////////////////////////////////////////////////////////////////
//
// File:     planar.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000111
// Version:  1.2
//
// Copyright (C) 1998-2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef PLANAR_H
#define PLANAR_H

#include "../../waveguide.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Planar
//
//   Uniform infinite planar layer, focusing on a single plane wave
//   from a non-interacting continuum.
//
//   Note: the convention for the sign of r_tm is such that
//   r_te = r_tm for theta=0.
//  
/////////////////////////////////////////////////////////////////////////////

class Planar : public MonoWaveguide
{ 
  public:

    Planar(Material& m);
  
    Complex eps_at(const Coord& coord) const
      {return core->eps();}
    
    Complex mu_at (const Coord& coord) const
      {return core->mu();}

    bool operator==(const Waveguide& w) const
      {return *core == *(dynamic_cast<const Planar*>(&w)->core);}
        
    static void set_kt(const Complex& kt_)
      {kt = kt_;}

    static Complex get_kt()
      {return kt;}
    
    void set_theta(Real theta_radians)
      {kt = (2*pi / global.lambda * core->n()) * sin(theta_radians);}

    Complex c1_size() const 
      {return 0.0;}

    void find_modes()
      {mode->kz = calc_kz(); mode->kz_bw = -mode->kz;}

  protected:

    // Transverse component of wavevector is same for all layers in stack
    // because of Snell's law.
    
    static Complex kt;

    Complex calc_kz();
};

    
    
#endif



