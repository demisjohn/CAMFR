
/////////////////////////////////////////////////////////////////////////////
//
// File:     mode.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990405
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef MODE_H
#define MODE_H

#include <string>
#include <iostream>
#include "defs.h"
#include "coord.h"
#include "field.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Mode
//
//   Mode represented by polarisation, kz and field profile.
//   Should be normalised in a suitable way.
//   Is not a pure virtual class, since it is used in the copy constructor
//   of waveguide.
//  
/////////////////////////////////////////////////////////////////////////////

class Mode
{
  public:

     Mode(const Polarisation pol_,
          const Complex& kz_,   const Complex& kz_bw_,
          const Complex A_=0.0, const Complex B_=0.0)
       : pol(pol_), kz(kz_), kz_bw(kz_bw_), A(A_), B(B_) {}

     virtual void normalise() {} // Should be called after constructor.
     
     virtual ~Mode() {}
     
     virtual Field field(const Coord& coord) const {Field f; return f;}
     // Note: field is calculated at z=0. Propgation is done in field.cpp.

     virtual Field bw_field(const Coord& coord) const {return field(coord);}
     // Tmp: will be corrected when dealing with non-reciproque media.
     virtual Complex get_kz()    const {return kz;}
     virtual Complex get_kz_bw() const {return kz_bw;}

     Complex n_eff()    const {return get_kz()   /2.0/pi*global.lambda;}
     Complex n_eff_bw() const {return get_kz_bw()/2.0/pi*global.lambda;}
          
     Complex E_cst() const {return A;}
     Complex H_cst() const {return B;}     
     
     bool operator==(const Mode& rhs) const;

     std::string repr() const;

     Polarisation pol;
     
  protected:

     Complex kz;
     Complex kz_bw;

     Complex A; // proportionality constant in Ez
     Complex B; // proportionality constant in Hz

     friend class Planar;
};

inline std::ostream& operator<<(std::ostream& s, const Mode& mode)
  {return s << mode.repr() << std::endl;}




/////////////////////////////////////////////////////////////////////////////
//
// Function objects for sorting modes.
//
/////////////////////////////////////////////////////////////////////////////

struct betasorter
{
    bool operator()(const Complex& beta_a, const Complex& beta_b)
    {
      return ( real(beta_a * beta_a) > real(beta_b * beta_b) );
    }
};

struct modesorter
{
    bool operator()(const Mode* a, const Mode* b)
    {
      if ( (a->pol == TE) && (b->pol != TE) )
        return true;

      if ( (a->pol == TM) && (b->pol != TM) )
        return false;

      const Complex kz_a = a->get_kz();
      const Complex kz_b = b->get_kz();
      
      return ( real(kz_a * kz_a) > real(kz_b * kz_b) );
    }
};



/////////////////////////////////////////////////////////////////////////////
//
// Function object for sorting Bloch modes according to their losses. 
//
/////////////////////////////////////////////////////////////////////////////

struct modesorter_bloch
{
    bool operator()(const Mode* a, const Mode* b)
    {
      const Complex kz_a = a->get_kz();
      const Complex kz_b = b->get_kz();

      if (abs(kz_a) < 1e-8)
        return false;

      if (abs(kz_b) < 1e-8)
        return true;

      if ((abs(imag(kz_a)) > 1e-6) || (abs(imag(kz_b)) > 1e-6))
      {
        if (abs(imag(kz_a)+imag(kz_b)) < 1e-6)
          return imag(kz_a) < imag(kz_b);
        
        return ( abs(imag(kz_a)) < abs(imag(kz_b)) );
      }

      return real(kz_a) > real(kz_b);
    }
};




/////////////////////////////////////////////////////////////////////////////
//
// Sorts the modes to agree with the order in BDM's software (for easier
// comparison):
//
//   - n  = 0: first all the TE modes, then the TM modes
//              normal order for radiation modes
//   - n != 0: reverse order for EH and HE radiation modes
//
/////////////////////////////////////////////////////////////////////////////

struct modesorter_BDM
{
    bool operator()(const Mode* a, const Mode* b)
    {
      if ( (a->pol == TE) && (b->pol != TE) )
        return true;

      if ( (a->pol == TM) && (b->pol != TM) )
        return false;

      const double ra2 = real(a->get_kz() * a->get_kz());
      const double rb2 = real(b->get_kz() * b->get_kz());
      
      if ( (ra2 < 0) && (rb2 < 0) ) {
        if ( (a->pol != TE) && (a->pol != TM) )
          return ( ra2 < rb2 );
        else
          return ( ra2 > rb2 );
      }

      if ( (ra2 > 0) && (rb2 > 0) )
        return ( ra2 > rb2 );

      if ( (ra2 > 0) && (rb2 < 0) )
        return true;
      else
        return false;
    }
};



#endif

