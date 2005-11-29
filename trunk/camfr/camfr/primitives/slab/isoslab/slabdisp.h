
/////////////////////////////////////////////////////////////////////////////
//
// File:     slabdisp.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010314
// Version:  1.2
// Copyright (C) 2000-2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SLABDISP_H
#define SLABDISP_H

#include "slab.h"
#include "../../../math/calculus/calculus.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: SlabDisp
//
//   Dispersion relation for modes of a stack of infinite planar layers
//   sandwiched between two metal plates parallel to the z-axis.
//
//   kt = kx in the region where k is kt_eps_mu.
//
/////////////////////////////////////////////////////////////////////////////

class SlabDisp : public Function1D<Complex>
{
  public:

    SlabDisp(const Expression& ex, const Complex& lambda,
             SlabWall* lowerwall=NULL, SlabWall* upperwall=NULL);

    SlabDisp(const std::vector<Material*>& materials,
             const std::vector<Complex>&   thicknesses,
             const Complex& lambda,
             SlabWall* lowerwall=NULL, SlabWall* upperwall=NULL);

    ~SlabDisp() {}
    
    Complex operator()(const Complex& kt);

    Complex get_kt_eps_mu() const {return kt_eps_mu;}

    std::vector<Complex> get_params() const;
    
    void set_params(const std::vector<Complex>& params);
    
  protected:

    std::vector<Complex> thicknesses;
    std::vector<Complex> eps;
    std::vector<Complex> mu;
    
    Complex lambda;

    Complex kt_eps_mu;
    
    SlabWall* lowerwall;
    SlabWall* upperwall;
};



#endif
