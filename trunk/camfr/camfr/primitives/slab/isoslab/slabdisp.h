
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
//   kt = kx in the region with minimum refractive index.
//
/////////////////////////////////////////////////////////////////////////////

class SlabDisp : public Function1D<Complex>
{
  public:

    SlabDisp(const Expression& ex,    Real lambda,
             SlabWall* leftwall=NULL, SlabWall* rightwall=NULL);

    SlabDisp(const std::vector<Material*>& materials,
             const std::vector<Complex>&   thicknesses,
             Real lambda, SlabWall* leftwall=NULL, SlabWall* rightwall=NULL);

    ~SlabDisp() {}
    
    Complex operator()(const Complex& kt);


    Complex get_min_eps_mu() const {return min_eps_mu;}

    std::vector<Complex> get_params() const;
    
    void set_params(const std::vector<Complex>& params);
    
  protected:

    std::vector<Complex> thicknesses;
    std::vector<Complex> eps;
    std::vector<Complex> mu;
    
    Real lambda;

    Complex min_eps_mu;

    SlabWall*  leftwall;
    SlabWall* rightwall;
};



#endif
