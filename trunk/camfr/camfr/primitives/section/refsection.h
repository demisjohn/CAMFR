
/////////////////////////////////////////////////////////////////////////////
//
// File:     refsection.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20021104
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef REFSECTION_H
#define REFSECTION_H

#include <vector>
#include "sectionmode.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: RefSection
//
//  Uniform reference section for use in the series expansion root finder.
//  Is not entirely equivalent to a Section1D with a UniformSlab, as here
//  we take different linear combinations of degenerate modes to get modes
//  which are nicely TE or TM with respect to the xy plane.
//  
/////////////////////////////////////////////////////////////////////////////

class RefSection : public SectionImpl
{
  public:

    RefSection(Material& mat, const Complex& width, const Complex& height) 
      : m(&mat), a(width), b(height) {uniform=true; core=m; M=1;}

    ~RefSection() {}
    
    bool operator== (const Waveguide& w) const
      {return this == &w;}
    
    std::vector<Material*> get_materials() const
      {std::vector<Material*> v; v.push_back(m); return v;}
    
    bool contains(const Material& mat) const
      {return m == &mat;}

    bool no_gain_present() const
      {return m->no_gain_present();}

    Complex eps_at(const Coord& c) const
      {return m->eps();}

    Complex get_eps() const
      {return m->eps();}

    Complex  mu_at(const Coord& c) const
      {return m->mu();}

    Complex get_mu() const
      {return m->mu();}

    Complex get_width() const
      {return a;}

    Complex get_height() const
      {return b;}

    void find_modes();

    void calc_overlap_matrices(Section2D* profile, 
                               cMatrix* O_EE, cMatrix* O_MM, 
                               cMatrix* O_EM, cMatrix* O_zz);

  protected:

    Material* m;
    Complex a;
    Complex b;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: RefSectionMode
//
/////////////////////////////////////////////////////////////////////////////

class RefSectionMode : public Mode
{
  public:

    RefSectionMode(Polarisation pol,  const Complex& kz,
                   const Complex& kx, const Complex& kx0,
                   const Complex& ky, const Complex& ky0, 
                   RefSection* geom);

    Field field(const Coord& coord) const;

  protected:

    Complex kx, kx0;
    Complex ky, ky0;
    RefSection* geom;
};



#endif
