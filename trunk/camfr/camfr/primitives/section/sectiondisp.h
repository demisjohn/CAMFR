
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectiondisp.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020125
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SECTIONDISP_H
#define SECTIONDISP_H

#include "../../stack.h"
#include "../slab/generalslab.h"
#include "../../math/calculus/calculus.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SectionDisp
//
//  Dispersion relation for a waveguide with a 2D cartesian cross section
//
/////////////////////////////////////////////////////////////////////////////

class SectionDisp : public ComplexFunction
{
  public:

    SectionDisp(Stack& _left, Stack& _right, Real _lambda);

    Complex operator()(const Complex& beta);

    vector<Complex> get_params() const;
    void set_params(const vector<Complex>&);
    
  protected:

    Stack* left;
    Stack* right;

    vector<Slab*> left_slabs;
    vector<Slab*> right_slabs;  

    vector<Complex*> left_d;
    vector<Complex*> right_d;
    
    Real lambda;
};



#endif
