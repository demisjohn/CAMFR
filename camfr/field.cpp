
/////////////////////////////////////////////////////////////////////////////
//
// File:     field.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990505
// Version:  1.1
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <iostream>
#include "waveguide.h"
#include "field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Field::operator+
//  
/////////////////////////////////////////////////////////////////////////////

Field Field::operator+ (const Field& b) const
{
  Field result;

  result.E1 = E1 + b.E1;  
  result.E2 = E2 + b.E2;
  result.Ez = Ez + b.Ez;

  result.H1 = H1 + b.H1;  
  result.H2 = H2 + b.H2;
  result.Hz = Hz + b.Hz;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Field::operator-
//  
/////////////////////////////////////////////////////////////////////////////

Field Field::operator- (const Field& b) const
{
  Field result;

  result.E1 = E1 - b.E1;  
  result.E2 = E2 - b.E2;
  result.Ez = Ez - b.Ez;

  result.H1 = H1 - b.H1;  
  result.H2 = H2 - b.H2;
  result.Hz = Hz - b.Hz;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Field::operator*
//  
/////////////////////////////////////////////////////////////////////////////

Field Field::operator* (const Complex& c) const
{ 
  Field result;

  result.E1 = E1 * c;  
  result.E2 = E2 * c;
  result.Ez = Ez * c;

  result.H1 = H1 * c;  
  result.H2 = H2 * c;
  result.Hz = Hz * c;

  return result;
}


/////////////////////////////////////////////////////////////////////////////
//
// Field::operator/
//  
/////////////////////////////////////////////////////////////////////////////

Field Field::operator/ (const Complex& c) const
{
  Field result;
  
  result.E1 = E1 / c;  
  result.E2 = E2 / c;
  result.Ez = Ez / c;

  result.H1 = H1 / c;  
  result.H2 = H2 / c;
  result.Hz = Hz / c;

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// Field::repr
//  
/////////////////////////////////////////////////////////////////////////////

std::string Field::repr() const
{
  std::ostringstream s;
  
  s << "E1="  << E1 << ", E2=" << E2 << ", Ez=" << Ez << std::endl;
  s << "H1="  << H1 << ", H2=" << H2 << ", Hz=" << Hz;
  
  return  s.str();
}



/////////////////////////////////////////////////////////////////////////////
//
// FieldExpansion::field
//  
/////////////////////////////////////////////////////////////////////////////

Field FieldExpansion::field(const Coord& coord) const
{ 
  Field f;
  
  for (unsigned int i=1; i<=wg->N(); i++)
  {
    Field f_i = wg->get_mode(i)->field(coord);
    Complex I_kz_d = I * wg->get_mode(i)->get_kz() * coord.z;
    
    f.E1 += fw(i) * f_i.E1 * exp(-I_kz_d)  +  bw(i) * f_i.E1 * exp(I_kz_d);
    f.E2 += fw(i) * f_i.E2 * exp(-I_kz_d)  +  bw(i) * f_i.E2 * exp(I_kz_d);
    f.Ez += fw(i) * f_i.Ez * exp(-I_kz_d)  -  bw(i) * f_i.Ez * exp(I_kz_d);

    f.H1 += fw(i) * f_i.H1 * exp(-I_kz_d)  -  bw(i) * f_i.H1 * exp(I_kz_d);
    f.H2 += fw(i) * f_i.H2 * exp(-I_kz_d)  -  bw(i) * f_i.H2 * exp(I_kz_d);
    f.Hz += fw(i) * f_i.Hz * exp(-I_kz_d)  +  bw(i) * f_i.Hz * exp(I_kz_d);
  }

  return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// FieldExpansion::propagate
//  
/////////////////////////////////////////////////////////////////////////////

FieldExpansion FieldExpansion::propagate(const Complex& z) const
{
  cVector fw_c(fw.rows(), fortranArray);
  cVector bw_c(bw.rows(), fortranArray);
  
  FieldExpansion f(wg, fw_c, bw_c);
  
  for (unsigned int i=1; i<=wg->N(); i++)
  {
    Complex I_kz_d = I * wg->get_mode(i)->get_kz() * z;

    fw_c(i) = fw(i) * exp(-I_kz_d);
    bw_c(i) = bw(i) * exp( I_kz_d);    
  }

  return f;
}



/////////////////////////////////////////////////////////////////////////////
//
// FieldExpansion::operator*
//  
/////////////////////////////////////////////////////////////////////////////

FieldExpansion FieldExpansion::operator* (const Complex& c) const
{
  cVector fw_c(fw); fw_c.makeUnique(); fw_c *= c;
  cVector bw_c(bw); bw_c.makeUnique(); bw_c *= c;
  
  return FieldExpansion(wg, fw_c, bw_c);
}

 

/////////////////////////////////////////////////////////////////////////////
//
// FieldExpansion::repr
//  
/////////////////////////////////////////////////////////////////////////////

std::string FieldExpansion::repr() const
{
  std::ostringstream s;
  
  s << "Core: " << wg->get_core()->n() << std::endl;
  s << "fw: " << fw << std::endl;
  s << "bw: " << bw << std::endl;

  return s.str();
}
