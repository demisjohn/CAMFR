
/////////////////////////////////////////////////////////////////////////////
//
// File:     cavity.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20000104
// Version:  1.0
//
// Copyright (C) 2000 Peter Bienstman - Ghent University
//
////////////////////////////////////////////////////////////////////////////

#ifndef CAVITY_H
#define CAVITY_H

#include "stack.h"
#include "math/calculus/function.h"

class Cavity; // forward declaration

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Sigma_lambda
//
//  A function object to calculate the singular value of a cavity as a
//  function of the wavelength.
//  
/////////////////////////////////////////////////////////////////////////////

class Sigma_lambda : public Function1D<Real>
{
  public:

    Sigma_lambda(Cavity* cavity_) : cavity(cavity_) {}

    Real operator()(const Real& lambda);

  protected:

    Cavity* cavity;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Sigma_n_imag
//
//  A function object to calculate the singular value of a cavity as a
//  function of the material gain.
//  
/////////////////////////////////////////////////////////////////////////////

class Sigma_n_imag : public Function1D<Real>
{
  public:

    Sigma_n_imag(Cavity* cavity_) : cavity(cavity_) {}

    Real operator()(const Real& lambda);

  protected:

    Cavity* cavity;
};



/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Cavity
//
//  A cavity sliced in half in an arbitrary plane to form a top and bottom
//  part.
//  
/////////////////////////////////////////////////////////////////////////////

class Cavity
{
  public:

    Cavity(Stack& bot, Stack& top);

    // Look for the sourceless modes in a given wavelength region.
    // If number=0, all modes are subsequently located using 'find_mode',
    // else only the highest 'number' modes are located.

    void find_modes_in_region(Real lambda_start, Real lambda_stop,
                              Real delta_lambda, unsigned int number=0,
                              Real n_imag_start=0.0, Real n_imag_stop=0.015,
                              unsigned int passes=1);

    // Locate the single sourceless mode confined in a given wavelength (um)
    // and gain region, using 'passes' number of passes.     

    void find_mode(Real lambda_start,     Real lambda_stop,
                   Real n_imag_start=0.0, Real n_imag_stop=0.015,
                   unsigned int passes=1);

    // Calculate sigma and dominant mode for sourceless problem with
    // current wavelength and gain. Also sets cavity field profile.

    Real calc_sigma(int* dominant_mode=NULL);

    // Calculates total field in the cavity after introduction of a 
    // general source.

    void set_source(const cVector& fw, const cVector& bw)
      {set_source_expansion(FieldExpansion(top->get_inc(),fw,bw));}

    // Calculates total field in the cavity after introduction of a 
    // current source.

    void set_source(const Coord& pos, const Coord& orientation)
      {set_source_expansion(dynamic_cast<MultiWaveguide*>(top->get_inc())
                  ->field_from_source(pos, orientation));}
    
    // Calculate field profile in cavity.
    
    Field field(const Coord& coord);

    // Get index profile.

    Complex eps_at(const Coord& c) const;
    Complex  mu_at(const Coord& c) const;
    Complex   n_at(const Coord& c) const {return sqrt(eps_at(c)/eps0);}

  protected:

    void set_source_expansion(const FieldExpansion& f);

    Stack* top;
    Stack* bot;

    Sigma_lambda sigma_lambda;
    Sigma_n_imag sigma_n_imag;

    Real current_sigma;
};



#endif






