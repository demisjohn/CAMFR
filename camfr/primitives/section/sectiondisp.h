
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
#include "../../math/calculus/calculus.h"
#include "../slab/generalslab.h"

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

    SectionDisp(Stack& _left, Stack& _right, Real _lambda, int _M, 
                bool symmetric = false);

    Complex operator()(const Complex& beta);

    vector<Complex> get_params() const;
    void set_params(const vector<Complex>&);
    
  protected:

    Complex calc_lapack (const Complex& beta);
    Complex calc_lapack2(const Complex& beta);
    Complex calc_arnoldi(const Complex& beta);

    Stack* left;
    Stack* right;

    bool symmetric;

    vector<Slab*> left_slabs;
    vector<Slab*> right_slabs;  

    vector<Complex*> left_d;
    vector<Complex*> right_d;
    
    Real lambda;

    int M;
};



#endif
