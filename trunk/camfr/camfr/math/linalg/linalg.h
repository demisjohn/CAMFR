
/////////////////////////////////////////////////////////////////////////////
//
// File:     linalg.h
// Authors:  Peter.Bienstman@rug.ac.be, Lieven.Vanholme@rug.ac.be
// Date:     20020207
// Version:  1.1
//
// Copyright (C) 1998-2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef LINALG_H
#define LINALG_H

#include "../../defs.h"
#include "blitz/array.h"

typedef blitz::Array<Complex,1> cVector;
typedef blitz::Array<Complex,2> cMatrix;
typedef blitz::Array<Complex,3> cHyperM;
typedef blitz::Array<Real,1>    rVector;
typedef blitz::Array<Real,2>    rMatrix;
typedef blitz::Array<Real,3>    rHyperM;
typedef blitz::Array<int,1>     iVector;
typedef blitz::Array<int,2>     iMatrix;
typedef blitz::Array<int,3>     iHyperM;

using blitz::fortranArray;

typedef enum{nrml, transp, herm} Op;

/////////////////////////////////////////////////////////////////////////////
//
// Front end for linear algebra manipulations. Uses Blitz C++ Array class
// and LAPACK Fortran numerical library.
//
// Using these routines assumes that the arrays are created using Fortran
// storage ordering. E.g. for a 10 by 10 complex matrix:
//
//    cMatrix(10,10,fortranArray) A;
//
// Take care that the matrix used to store the results of an operation
// has the correct dimensions.
//
// Currently, only operations on complex matrices are provided. For real
// matrices, 'Complex' should be changed to 'Real' in C++ and the Fortran
// LAPACK prefix should change from 'Z' to 'D'.
//
// For performance reasons, write e.g.
//   C.reference(multiply(A,B))
// instead of
//   C = multiply(A,B)
// This avoids unneccesary copying of the matrix data.
//  
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// returns transpose or hermitian conjugate of A
//
/////////////////////////////////////////////////////////////////////////////

cMatrix transpose(const cMatrix& A);
cMatrix herm_conj(const cMatrix& A);



/////////////////////////////////////////////////////////////////////////////
//
// Replaces a square matrix A by its transpose or hermitian conjugate.
//
/////////////////////////////////////////////////////////////////////////////

void transpose_self(cMatrix* A);
void herm_conj_self(cMatrix* A);



/////////////////////////////////////////////////////////////////////////////
//
// Matrix with vector multiplication: multiply(A,x).
//  
/////////////////////////////////////////////////////////////////////////////

cVector multiply(const cMatrix& A, const cVector& x, Op a=nrml);



/////////////////////////////////////////////////////////////////////////////
//
// Matrix multiplications: multiply(A,B)
//                         multiply(Op(A),Op(B))
//                         multiply(A,B,C)
//                         multiply(A,B,C,D)
//
// Note: don't write A*B for matrix multiplications, this will do
// element-wise multiplication !
//  
/////////////////////////////////////////////////////////////////////////////

cMatrix multiply(const cMatrix& A, const cMatrix& B, Op a=nrml, Op b=nrml);

inline cMatrix multiply(const cMatrix& A, const cMatrix& B, const cMatrix& C)
    {return multiply(A,multiply(B,C));}

inline cMatrix multiply(const cMatrix& A, const cMatrix& B,
                        const cMatrix& C, const cMatrix& D)
    {return multiply(A,multiply(B,C,D));}

inline cVector multiply(const cMatrix& A, const cMatrix& B, const cVector& x)
    {return multiply(A,multiply(B,x));}



/////////////////////////////////////////////////////////////////////////////
//
// Returns solution X of system of linear equations A*X=B.
//
//   B can contain more than one column.  
//   Note: don't use explicit matrix inversion here, since that is slower.
//
//   solve_x uses expert Lapack driver to improve stability
//   solve_svd solves ipoorly condtioned system using SVD.
//
//   The 'sym' variants expect a symmetrix matrix A.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix       solve(const cMatrix& A, const cMatrix& B);
cMatrix     solve_x(const cMatrix& A, const cMatrix& B);
cMatrix   solve_svd(const cMatrix& A, const cMatrix& B);
cMatrix   solve_sym(const cMatrix& A, const cMatrix& B);
cMatrix solve_sym_x(const cMatrix& A, const cMatrix& B);



/////////////////////////////////////////////////////////////////////////////
//
// Computes eigenvalues and/or eigenvectors of matrix A.
//
//   The i-th column of matrix 'eigenvectors' corresponds to the i-th
//   eigenvalue.
//
//   eigenvalues_x uses expert Lapack driver to improve stability.
//
/////////////////////////////////////////////////////////////////////////////

cVector eigenvalues  (const cMatrix& A, cMatrix* eigenvectors=NULL);
cVector eigenvalues_x(const cMatrix& A, cMatrix* eigenvectors=NULL);



/////////////////////////////////////////////////////////////////////////////
//
// Computes eigenvalues and/or eigenvectors of the generalized 
// eigenproblem Ax = lambda Bx, where lambda = alpha / beta.
//
/////////////////////////////////////////////////////////////////////////////

void gen_eigenvalues(const cMatrix& A, const cMatrix& B,
                     cVector* alpha, cVector* beta,
                     cMatrix* eigenvectors=NULL);



/////////////////////////////////////////////////////////////////////////////
//
// Computes the SVD decomposition U.Sigma.Vh of matrix A.
//
/////////////////////////////////////////////////////////////////////////////

rVector svd(const cMatrix& A, cMatrix* Vh=NULL, cMatrix* U=NULL);


     
/////////////////////////////////////////////////////////////////////////////
//
// Inverts a matrix.
//
//   Note: mostly not needed explicitly. Use routines to solve a linear
//   system instead, since that is faster.
//
//   invert_x uses expert Lapack driver to improve stability.
//   invert_svd inverts a poorly condtioned matrix using SVD.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix invert    (const cMatrix& A);
cMatrix invert_x  (const cMatrix& A);
cMatrix invert_svd(const cMatrix& A);



/////////////////////////////////////////////////////////////////////////////
//
// Determinant.
//
//   See LAPACK user guide for details on banded storage.
//
/////////////////////////////////////////////////////////////////////////////

Complex determinant     (const cMatrix& A);
Complex determinant_band(const cMatrix& A, int rows, int kl, int ku);



/////////////////////////////////////////////////////////////////////////////
//
// Do LU decomposition of matrix A.
//
//   Low level routine. Mostly not needed to use explicitly.
//
/////////////////////////////////////////////////////////////////////////////

void LU(const cMatrix& A, cMatrix* LU, iVector* P);



/////////////////////////////////////////////////////////////////////////////
//
// Solve linear system op(A)*X=B starting from LU decomposition of A.
//
//   op = nrml, transp or herm
//
//   Low level routine. Mostly not needed to use explicitly, since the
//   routine 'solve' has the same functionality for op = nrml.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix LU_solve(const cMatrix& LU, const iVector& P,
                 const cMatrix& B, Op op=nrml);



#endif



