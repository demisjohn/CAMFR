
/////////////////////////////////////////////////////////////////////////////
//
// File:     defs.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19990201
// Version:  1.01
//
// Copyright (C) 1998,1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef DEFS_H
#define DEFS_H

#include <math.h>
#include <complex>
#include <string>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//
// General conventions
//
//  - Time dependence is given by exp(jwt).
//  - All array and vector indexing starts at 1.
//  - Arrays are laid out in Fortran storage order.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// typedefs
//
/////////////////////////////////////////////////////////////////////////////

typedef double Real;
typedef complex<Real> Complex;
typedef enum {Plus, Min} Limit;
typedef enum {ADR, track} Solver;
typedef enum {normal, extra, SVD} Stability;
typedef enum {T_T, S_T, S_S} Field_calc;  
typedef enum {GEV, T} Bloch_calc;



/////////////////////////////////////////////////////////////////////////////
//
// typedefs for polarisation
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {unknown, TEM, TE, TM, HE, EH, TE_TM} Polarisation;
extern const string Pol_string[11];
ostream& operator<< (ostream& s, const Polarisation& pol);



/////////////////////////////////////////////////////////////////////////////
//
// constants
//
/////////////////////////////////////////////////////////////////////////////

const Real pi   = 3.141592653589793238462643383279502884197;
const Real c    = 2.997925e8;
const Real mu0  = 4*pi*1e-7;
const Real eps0 = 1.0/mu0/c/c;

const Complex I(0.0, 1.0);

Real machine_eps();



/////////////////////////////////////////////////////////////////////////////
//
// STRUCT: Global
//
//   Groups global variables and numerical parameters
//
/////////////////////////////////////////////////////////////////////////////

class Material; // forward declaration, see material.h
 
struct Global
{
    // Wavelength.
    Real lambda;

    // Number of modes retained in eigenmode expansion. 
    unsigned int N;

    // Polarisation, if applicable for geometries with uncouped TE/TM modes.
    Polarisation polarisation;

    // The material whose gain will be adjusted to get amplitude resonance.
    Material* gain_mat;

    // Complex solver to use: ADR or root tracking.
    Solver solver;
    
    // Stability mode: normal, extra matrix equilibration or SVD.
    Stability stability;

    // Precision factor for resolution in n_eff in finding guided modes
    // on the coordinate axes. Higher number is more precise.
    unsigned int precision;

    // Precision enhancement factor for finding guided modes. If larger
    // than 1, a search with precision of 'precision*precision_enhancement'
    // will be done around the guided modes found with precision of
    // 'precision'. Usefull for modes that are nearly degenerate.
    unsigned int precision_enhancement;

    // Specifies the width of the region that has to be searched around the
    // original guided mode estimates, if precision_enhancement is larger
    // than one. This region is bounded by (1-dx_enhanced)*x_estimate and
    // (1+dx_enhanced)*x_estimate.
    Real dx_enhanced;

    // Precision factor for resolution in n_eff in finding radiation modes
    // on the coordinate axes. Higher number is more precise.
    unsigned int precision_rad;

    // Precision factor for number of steps in complex countour root solver.
    // Higher number is more precise.
    unsigned int C_steps;

    // Relative factor in determining the upper right of the rectangle in
    // the kt plane that is searched for complex zeros.
    // If e.g. chosen to be (2,3), upperright will be
    // (2*real_upperright_default, 3*imag_upperright_default)
    Complex C_upperright;

    // When changing the geometry parameters or wavelength, the modes are
    // found starting from those of the previous structure.
    bool sweep_from_previous;

    // Initial number of steps taken in tracing root solver.
    // higher number is more precise.
    unsigned int sweep_steps;

    // Coarse intermediate precision for root finding during tracing root
    // solver.
    Real eps_trace_coarse;

    // Determines if during root tracing all roots are traced together or
    // in chunks.
    bool chunk_tracing;

    // If the amplitude of an increasing exponential is smaller than this
    // number, the result is put to zero.
    Real unstable_exp_threshold;

    // Determines if the field profiles in stacks are calculated with S or
    // T type excitations.
    Field_calc field_calc;

    // Determines which algorithm is used to calculate Bloch modes.
    Bloch_calc bloch_calc;

    // TMP switch that tells CAMFR is the local modes are orthogonal.
    // CAMFR should later figure this out by itself.
    bool orthogonal;
};

extern Global global;



/////////////////////////////////////////////////////////////////////////////
//
// out_of_memory
//
//   error handler
//
/////////////////////////////////////////////////////////////////////////////

void out_of_memory();



#endif



