
/////////////////////////////////////////////////////////////////////////////
//
// File:     lintest.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19980524
// Version:  1.0
//
// Copyright (C) 1998 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// main function for testing some linear algebra functions
//
/////////////////////////////////////////////////////////////////////////////

int main()
{  
  //
  // multiplications
  //

  cMatrix A(3,2,fortranArray);

  A(1,1) = 3; A(1,2) = 4;
  A(2,1) = 1; A(2,2) = 3;
  A(3,1) = 1; A(3,2) = 3;
  
  cout << "A: " << A << endl;

  cVector x(2,fortranArray);

  x(1) = 1;
  x(2) = 2;

  cout << "x: " << x << endl;

  cVector y(3,fortranArray);
  y = multiply(A,x);

  cout << "y=A*x: " << y << endl << endl;
  
  cMatrix B(2,3,fortranArray);

  B(1,1) = 1; B(1,2) = 2; B(1,3) = 5;
  B(2,1) = 1; B(2,2) = 3; B(2,3) = 6;
  
  cout << "B: " << B << endl;

  cMatrix C(3,3,fortranArray);
  C = multiply(A,B);

  cout << "C=A*B: " << C << endl;

  //
  // linear systems
  //
  
  cMatrix D(3,3,fortranArray);

  D(1,1) = 1; D(1,2) = 2; D(1,3) = 5;
  D(2,1) = 6; D(2,2) = 3; D(2,3) = 6;
  D(3,1) = 1; D(3,2) = 3; D(3,3) = 2;

  cout << "D: " << D << endl;;
  
  cout << "solve D*X=A: " << solve(D,A) << endl;

  cout << "solve again (expert): "<< solve_x(D,A) << endl;

  cout << "low level:" << endl;

  cMatrix lu(3,3,fortranArray);
  iVector p(3,fortranArray);

  LU(D,&lu,&p);

  cout << LU_solve(lu,p,A,transp) << endl;

  //
  // eigenvalues
  //
  
  cMatrix eigenvectors(3,3,fortranArray);
  
  cout << "eigenvalues of D: " << endl;
  cout << eigenvalues(D,&eigenvectors) << endl << endl;
  cout << "eigenvectors of D: " << eigenvectors << endl;

  cout << "eigenvalues of D (expert): " << endl;
  cout << eigenvalues_x(D,&eigenvectors) << endl << endl;
  cout << "eigenvectors of D (expert): " << eigenvectors << endl;

  cout << "eigenvalues of D (again): " << endl;
  cout << eigenvalues(D) << endl;

  //
  // SVD
  //

  cMatrix Vh(B.columns(),B.columns(),fortranArray);
  cMatrix  U(B.rows(),   B.rows(),   fortranArray);

  cout << endl << "SVD of B: " << endl;
  cout << svd(B,&Vh,&U) << endl << endl;
  cout << "U : " << U << endl;
  cout << "Vh : " << Vh << endl;

  //
  // Inverse
  //

  cout << "Inv(D) : " << invert(D) << endl;

  cout << "Inv_svd(D) : " << invert_svd(D) << endl;
  
  return 0;  
}

  
