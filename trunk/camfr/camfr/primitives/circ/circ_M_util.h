
/////////////////////////////////////////////////////////////////////////////
//
// File:     circ_M_util.h
// Author:   Mihai Ibanescu (michel@mit.edu)
// Date:     20020402
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CIRC_M_UTIL_H
#define CIRC_M_UTIL_H

#include "../../defs.h"
#include "../../material.h"
#include "../../mode.h"
#include "../../math/calculus/calculus.h"
#include "../../math/bessel/bessel.h"
#include "../../math/linalg/linalg.h"

using std::vector;

#include "../../util/vectorutil.h"

/////////////////////////////////////////////////////////////////////////////
//
// Functions for calculating transfer matrices (for the transfer matrix method)
// and for calculating the fields a given point once we know the vector of amplitudes
//
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//
// Calculates transfer matrix for interface between mat1 and mat2 at a radius r
// Uses kz as independent variable (for the moment).
//
/////////////////////////////////////////////////////////////////////////////

cMatrix transfer_matrix(const Complex& r, 
			const Material& mat1, 
			const Material& mat2,
			const Complex& kz,
			const int& circ_order,
			const bool& scaling = false);

/////////////////////////////////////////////////////////////////////////////
//
// Calculates field matrix 
// To get the actual field one multiplies the field matrix with the amplitude vector
// The 6-component vector of field components has the ordering [Ez Ephi Er Hz Hphi Hr]
//
/////////////////////////////////////////////////////////////////////////////

cMatrix field_matrix(const Complex& r, 
		     const Material& mat,
		     const Complex& kz,
		     const int& circ_order,
		     const bool& scaling = false);


#endif
