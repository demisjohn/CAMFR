
/////////////////////////////////////////////////////////////////////////////
//
// File:     field.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990505
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FIELD_H
#define FIELD_H

#include <string>
#include "defs.h"
#include "coord.h"
#include "math/linalg/linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Field
//
//   Field values at a fixed point in space, in rectangular or circular
//   coordinates.
//  
/////////////////////////////////////////////////////////////////////////////

class Field
{
  public:

    Field() : E1(0.0), E2(0.0), Ez(0.0), H1(0.0), H2(0.0), Hz(0.0) {}
  
    Complex E1; // Ex or Er
    Complex E2; // Ey or Ephi
    Complex Ez;

    Complex H1; // Hx or Hr
    Complex H2; // Hy or Hphi
    Complex Hz;

    Complex S1() const {return E2*conj(Hz) - Ez*conj(H2);}
    Complex S2() const {return Ez*conj(H1) - E1*conj(Hz);}
    Complex Sz() const {return E1*conj(H2) - E2*conj(H1);} 

    Real abs_E() const
      {return sqrt(real(E1*conj(E1) + E2*conj(E2) + Ez*conj(Ez)));}

    Real abs_H() const
      {return sqrt(real(H1*conj(H1) + H2*conj(H2) + Hz*conj(Hz)));}

    Real abs_S() const
      {return sqrt(real(S1()*conj(S1()) + S2()*conj(S2()) + Sz()*conj(Sz())));}

    Field operator+ (const Field&)   const;
    Field operator- (const Field&)   const;
    Field operator* (const Complex&) const;
    Field operator/ (const Complex&) const;

    Field& operator+= (const Field& f)   {*this = *this+f; return *this;}
    Field& operator-= (const Field& f)   {*this = *this-f; return *this;}
    Field& operator*= (const Complex& c) {*this = *this*c; return *this;}
    Field& operator/= (const Complex& c) {*this = *this/c; return *this;}

    std::string repr() const;
};

inline std::ostream& operator<<(std::ostream& s, const Field& f)
  {return s << f.repr() << std::endl;}

inline Field operator*(const Complex& c, const Field& f) {return f*c;}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: FieldExpansion
//
//   Description of a field as coefficients of eigenmode expansion.
//  
/////////////////////////////////////////////////////////////////////////////

class Waveguide; // forward declaration - see waveguide.h

class FieldExpansion
{
  public:

    FieldExpansion(Waveguide* wg_, const cVector& fw_, const cVector& bw_);
    FieldExpansion(const FieldExpansion& f);

    Field field(const Coord& coord) const;

    FieldExpansion propagate(const Complex& z) const;

    Waveguide* wg;

    cVector fw; //  forward propagating field
    cVector bw; // backward propagating field

    FieldExpansion operator* (const Complex&) const;

    std::string repr() const;
};

inline std::ostream& operator<<(std::ostream& s, const FieldExpansion& f)
  {return s << f.repr() << std::endl;}

inline FieldExpansion operator*(const Complex& c, const FieldExpansion& f)
  {return f*c;}



#endif
