
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabwall.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000927
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SLABWALL_H
#define SLABWALL_H

#include "../../../stack.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabWall
//
//   Transverse boundary condition for Slabs.
//   Contains a function to create a starting field satisfying the boundary
//   condition, as well as a function to test the error of a field when
//   fullfilling the condition. 
//  
////////////////////////////////////////////////////////////////////////////

class SlabWall
{
  public:

    SlabWall() {}
    virtual ~SlabWall() {}

    virtual Complex get_R12() const=0;
    virtual void get_start_field(Complex* in, Complex* out) const=0;
    virtual Complex get_error(const Complex& in, const Complex& out) const=0;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabWallMixed
//
//   Transverse boundary condition for Slabs a.in_wave + b.out_wave = 0.
//   Reflection coefficient is r = -a/b.
//   Defaults to electric (r = -1) wall.
//  
////////////////////////////////////////////////////////////////////////////

class SlabWallMixed : public SlabWall
{
  public:

    SlabWallMixed(const Complex& a_=1.0, const Complex& b_=1.0)
      : a(a_), b(b_) {}

    Complex get_R12() const
      {return -a/b;}
    
    void get_start_field(Complex* in, Complex* out) const
      {*in=b; *out=-a;}   
        
    Complex get_error(const Complex& in, const Complex& out) const
      {return a*in + b*out;}
     
  protected:

    Complex a, b;
};

extern const SlabWallMixed slab_E_wall;
extern const SlabWallMixed slab_H_wall;
extern const SlabWallMixed slab_open_wall;



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabWall_TBC
//
//   Transparent boundary condition. Reflectionless for a given angle, i.e.
//   for a given wavevector component kx_0 perpendicular to the wall.
//  
////////////////////////////////////////////////////////////////////////////

class SlabWall_TBC : public SlabWall
{
  public:

    SlabWall_TBC(const Complex& kx_0_, const Material& m_)
      : kx_0(kx_0_), m(&m_) {}

    Complex get_R12() const;
    void get_start_field(Complex* in, Complex* out) const;    
    Complex get_error(const Complex& in, const Complex& out) const;
     
  protected:

    const Complex kx_0;
    const Material* m;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabWall_PC
//
//   Infinite number of periods of a stack as a Photonic Crystal
//   transverse boundary condition.
//   Condition is of course dependent on the current value of Planar::kt.
//  
////////////////////////////////////////////////////////////////////////////

class SlabWall_PC : public SlabWall
{
  public:

    SlabWall_PC(const Expression& e) : s(e)
      {
        if (s.get_inc() != s.get_ext())
        {
          std::cerr << "Currently, incidence and exit media should match in "
                    << "SlabWall_PC." << std::endl;
          exit (-1);
        }

        if (!e.no_gain_present())
          std::cout << "Warning: results might be incorrect for gain in "
                    << "SlabWall_PC." << std::endl;
      }  

    Complex get_R12() const;
    
    void get_start_field(Complex* in, Complex* out) const
      {*in=1; *out=get_R12();}
    
    Complex get_error(const Complex& in, const Complex& out) const
      {return -get_R12()*in + out;}
    
  protected:

    mutable MonoStack s;
};



#endif
