
/////////////////////////////////////////////////////////////////////////////
//
// File:     coord.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991104
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "sstream"
#include "coord.h"

/////////////////////////////////////////////////////////////////////////////
//
// Coord::operator+
//  
/////////////////////////////////////////////////////////////////////////////

Coord Coord::operator+ (const Coord& b) const
{
  Coord result;

  result.c1 = c1 + b.c1;  
  result.c2 = c2 + b.c2;
  result.z  = z  + b.z;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Coord::operator-
//  
/////////////////////////////////////////////////////////////////////////////

Coord Coord::operator- (const Coord& b) const
{
  Coord result;

  result.c1 = c1 - b.c1;  
  result.c2 = c2 - b.c2;
  result.z  = z  - b.z;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Coord::operator*
//  
/////////////////////////////////////////////////////////////////////////////

Coord Coord::operator* (const Complex& c) const
{ 
  Coord result;

  result.c1 = c * c1;  
  result.c2 = c * c2;
  result.z  = c * z;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Coord::operator/
//  
/////////////////////////////////////////////////////////////////////////////

Coord Coord::operator/ (const Complex& c) const
{ 
  Coord result;

  result.c1 = c1 / c;  
  result.c2 = c2 / c;
  result.z  = z  / c;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Coord::operator<
//  
/////////////////////////////////////////////////////////////////////////////

bool Coord::operator< (const Coord& c) const
{
  if (abs(c1) < abs(c.c1))
    return true;
  if (abs(c1) > abs(c.c1))
    return false;
  else if ( (c1_limit==Min) && (c.c1_limit==Plus) )
    return true;
  

  if (abs(c2) < abs(c.c2))
    return true;
  if (abs(c2) > abs(c.c2))
    return false;
  else if ( (c2_limit==Min) && (c.c2_limit==Plus) )
    return true;

  if (abs(z) < abs(c.z))
    return true;
  if (abs(z) > abs(c.z))
    return false;
  else if ( (z_limit==Min) && (c.z_limit==Plus) )
    return true;

  return false;
}



/////////////////////////////////////////////////////////////////////////////
//
// Coord::repr
//  
/////////////////////////////////////////////////////////////////////////////

string limit(Limit l) {return (l==Plus) ? "+" : "-";}

string Coord::repr() const
{
  ostringstream s;

  s << "(" << c1 << limit(c1_limit) << ","
           << c2 << limit(c2_limit) << ","
           <<  z << limit( z_limit) << ")";

  return s.str();
}
