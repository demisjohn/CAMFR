
/////////////////////////////////////////////////////////////////////////////
//
// File:     coord.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991104
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef COORD_H
#define COORD_H

#include <string>
#include "defs.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Coord
//
//   Coord of a point in space, in rectangular or circular coordinates
//   support for left or right limits.
//
/////////////////////////////////////////////////////////////////////////////

class Coord
{
  public:

    Coord(const Complex& i_c1       = 0.0,
          const Complex& i_c2       = 0.0,
          const Complex& i_z        = 0.0,
                Limit    i_c1_limit = Min,
                Limit    i_c2_limit = Min,
                Limit    i_z_limit  = Min)
      : c1(i_c1),             c2(i_c2),             z(i_z),
        c1_limit(i_c1_limit), c2_limit(i_c2_limit), z_limit(i_z_limit) {};
    
    Complex c1; // x or r 
    Complex c2; // y or phi
    Complex z;

    Limit c1_limit;
    Limit c2_limit;
    Limit z_limit;
    
    Coord operator+ (const Coord&)   const;
    Coord operator- (const Coord&)   const;
    Coord operator* (const Complex&) const;
    Coord operator/ (const Complex&) const;

    Coord& operator+= (const Coord& co)  {*this = *this+co; return *this;}
    Coord& operator-= (const Coord& co)  {*this = *this-co; return *this;}
    Coord& operator*= (const Complex& c) {*this = *this*c;  return *this;}
    Coord& operator/= (const Complex& c) {*this = *this/c;  return *this;}

    bool operator< (const Coord&) const;

    string repr() const; 
};

inline Coord operator*(const Complex& c, const Coord& co) {return co*c;}

inline ostream& operator<<(ostream& s, const Coord& c)
  {return s << c.repr() << endl;}



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: CoordStretcher
//
//   Performs a uniform strecthing of a real coord to complex coord.
//  
/////////////////////////////////////////////////////////////////////////////

class CoordStretcher
{
  public:

    CoordStretcher(Real c_start_=0, const Complex& alpha_=1)
      : c_start(c_start_), alpha(alpha_) {};
    
    Complex operator()(Real c)
      {return (c < c_start) ? c : c_start + alpha*(c-c_start);}

  protected:

    const Real    c_start;
    const Complex alpha;    
};



/////////////////////////////////////////////////////////////////////////////
//
// Function object for sorting parts of coordinates using their real parts.
//
/////////////////////////////////////////////////////////////////////////////

struct RealSorter
{
    bool operator()(const Complex& a, const Complex& b)
    {
      return ( real(a) < real(b) );
    }
};



#endif
