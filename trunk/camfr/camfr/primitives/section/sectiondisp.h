
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

    SectionDisp(Stack& _left, Stack& _right, const Complex& _lambda, int _M, 
                bool symmetric = false);

    Complex operator()(const Complex& kt2);

    std::vector<Complex> get_params() const;
    void set_params(const std::vector<Complex>&);
    
  protected:

    Complex calc_split();
    Complex calc_global();
    Complex calc_T();

    Stack* left;
    Stack* right;
    Stack* stack;

    bool symmetric;
    bool metallic;

    Complex kt_eps_mu;

    std::vector<Slab*> left_slabs;
    std::vector<Slab*> right_slabs;  
    std::vector<Slab*> slabs;  

    std::vector<Complex*> left_d;
    std::vector<Complex*> right_d;
    std::vector<Complex>  d;
    
    Complex lambda;

    int M;
};



#endif
