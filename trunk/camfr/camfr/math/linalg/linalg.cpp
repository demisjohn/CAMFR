
/////////////////////////////////////////////////////////////////////////////
//
// File:     linalg.cpp
// Authors:  Peter.Bienstman@rug.ac.be, Lieven.Vanholme@rug.ac.be
// Date:     20020207
// Version:  1.1
//
// Copyright (C) 1998-2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "linalg.h"

/////////////////////////////////////////////////////////////////////////////
//
// Returns transpose or hermitian conjugate of A.
//
/////////////////////////////////////////////////////////////////////////////

cMatrix transpose(const cMatrix& A)
{
  cMatrix A_trans(A.columns(), A.rows(), fortranArray);

  for (int i=1; i<=A.rows(); i++)
    for (int j=1; j<=A.columns(); j++)
      A_trans(j,i) = A(i,j);

  return A_trans;
}

cMatrix herm_conj(const cMatrix& A)
{
  cMatrix A_herm_conj(A.columns(), A.rows(), fortranArray);

  for (int i=1; i<=A.rows(); i++)
    for (int j=1; j<=A.columns(); j++)
      A_herm_conj(j,i) = conj(A(i,j));

  return A_herm_conj;
}



/////////////////////////////////////////////////////////////////////////////
//
// Replaces a square matrix A by its transpose or hermitian conjugate.
//
/////////////////////////////////////////////////////////////////////////////

void transpose_self(cMatrix* A)
{
  const int N = A->rows();
  if (N != A->cols())
  {
    py_error("Error: transpose_self only implemented for square matrices.");
    exit (-1);
  }

  Complex tmp;
  for (int i=1; i<=N; i++)
    for (int j=i+1; j<=N; j++)
    {
      tmp = (*A)(i,j);
      (*A)(i,j) = (*A)(j,i);
      (*A)(j,i) = tmp;
    }
}

void herm_conj_self(cMatrix* A)
{
  const int N = A->rows();
  if (N != A->cols())
  {
    py_error("Error: herm_conj_self only implemented for square matrices.");
    exit (-1);
  }
  
  Complex tmp;
  for (int i=1; i<=N; i++)
  {
    (*A)(i,i) = conj((*A)(i,i));
    for (int j=i+1; j<=N; j++)
    {
      tmp = conj((*A)(i,j));
      (*A)(i,j) = conj((*A)(j,i));
      (*A)(j,i) = tmp;
    }
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// multiply(A,x)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgemv)
  (const char*,const int&,const int&,const Complex&,const Complex*,
   const int&,const Complex*,const int&,const Complex&,const Complex*,
   const int&);

cVector multiply(const cMatrix& A, const cVector& x, Op a)
{
  // Set dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int x_rows = x.rows();

  if (    ( (a == nrml) && (A_cols != x_rows) )
       || ( (a != nrml) && (A_rows != x_rows) )  )
  {
    std::ostringstream s;
    s << "Error: dimensions for matrix multiplication don't match : ";
    s << "[" << A_rows << "," << A_cols << "]*[" << x_rows << "].";
    py_error(s.str());
    exit (-1);
  }

  // Prepare op strings.

  char op_a[] = "N";
  if (a == transp)
    op_a[0] = 'T';
  if (a == herm)
    op_a[0] = 'C';

  // Calculate product.
  
  cVector y(A_rows,fortranArray);
  
  F77NAME(zgemv)(op_a,A_rows,A_cols,1,A.data(),A_rows,
                 x.data(),1,0,y.data(),1);
  
  return y;
}



/////////////////////////////////////////////////////////////////////////////
//
// multiply(A,B)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgemm)
  (const char*,const char*,const int&,const int&,const int&,const Complex&,
   const Complex*,const int&,const Complex*,const int&,const Complex&,
   const Complex*,const int&);

cMatrix multiply(const cMatrix& A, const cMatrix& B, Op a, Op b)
{
  // Set dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (    ( (a == nrml) && (b == nrml) && (A_cols != B_rows) )
       || ( (a != nrml) && (b == nrml) && (A_rows != B_rows) )
       || ( (a == nrml) && (b != nrml) && (A_cols != B_cols) )
       || ( (a != nrml) && (b != nrml) && (A_rows != B_cols) )   )
  {
    std::ostringstream s;
    s << "Error: dimensions for matrix multiplication don't match : ";
    s << "[" << A_rows << "," << A_cols << "]*["
      << B_rows << "m" << B_cols << "].";
    py_error(s.str());
    exit (-1);
  }

  // Prepare op strings.

  char op_a[] = "N";
  if (a == transp)
    op_a[0] = 'T';
  if (a == herm)
    op_a[0] = 'C';

  char op_b[] = "N";
  if (b == transp)
    op_b[0] = 'T';
  if (b == herm)
    op_b[0] = 'C';

  // Calculate product.
  
  cMatrix C(A_rows,B_cols,fortranArray);

  // Blitz implementation for op's=nrml:
  //
  // firstIndex  i;
  // secondIndex j;
  // thirdIndex  k;
  //
  // C = sum(A(i,k) * B(k,j), k);

  // BLAS Fortran routine:
  
  F77NAME(zgemm)(op_a,op_b,A_rows,B_cols,A_cols,1,A.data(),A_rows,
                 B.data(),B_rows,0,C.data(),A_rows);
  
  return C;
}



/////////////////////////////////////////////////////////////////////////////
//
// solve(A,B)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgesv)
  (const int&,const int&,const Complex*,const int&,const int*,
   const Complex*,const int&,int&);

cMatrix solve(const cMatrix& A, const cMatrix& B)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  // Make copy of matrices, since LAPACK overwrites them.

  cMatrix A_LU(A_rows,A_cols,fortranArray);
  A_LU = A;
  
  cMatrix B_X(B_rows,B_cols,fortranArray);
  B_X = B;
  
  iVector P(A_rows,fortranArray); // pivot indices

  // Solve system.
  
  int info;
  
  F77NAME(zgesv)(A_rows,B_cols,A_LU.data(),A_rows,
                 P.data(),B_X.data(),B_rows,info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }
  
  return B_X;
}



/////////////////////////////////////////////////////////////////////////////
//
// solve_x(A,B)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgesvx)
  (const char*,const char*,const int&,const int&,const Complex*,const int&,
   const Complex*,const int&,const int*,const char*,const Real*,const Real*,
   const Complex*,const int&,const Complex*,const int&,Real&,const Real*,
   const Real*,const Complex*,const Real*,int&);

cMatrix solve_x(const cMatrix& A, const cMatrix& B)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  // Make copy of matrices, since LAPACK overwrites them.

  cMatrix A_bis(A_rows,A_cols,fortranArray);
  A_bis = A;
  
  cMatrix B_bis(B_rows,B_cols,fortranArray);
  B_bis = B;

  // Prepare output arguments.
  
  cMatrix AF(A_rows,A_cols,fortranArray);  // factored A
  iVector P(A_rows,fortranArray);          // pivot indices
  char equi[1];                            // equilibration performed?
  rVector R(A_rows,fortranArray);          // row scaling factors
  rVector C(A_cols,fortranArray);          // column scaling factors
  cMatrix X(B_rows,B_cols,fortranArray);   // solution
  Real cond;                               // condition number
  rVector ferr(B_rows,fortranArray);       // forward error bound
  rVector berr(B_rows,fortranArray);       // backward error bound

  // Create workspace.
  
  cVector work (2*A_rows,fortranArray);
  rVector work2(2*A_rows,fortranArray);

  // Solve system.
  
  int info;

  F77NAME(zgesvx)("E","N",A_rows,B_cols,A_bis.data(),A_rows,AF.data(),A_rows,
                  P.data(),equi,R.data(),C.data(),B_bis.data(),B_rows,
                  X.data(),B_rows,cond,ferr.data(),berr.data(),work.data(),
                  work2.data(),info);

  // Diagnostic output:

  //cout << "Equi: " << equi[0] << endl;
  //cout << "Cond: " << cond << endl;
  //cout << "Ferr: " << ferr << endl;
  //cout << "Berr: " << berr << endl;

  if (info == 0)
    return X;
  
  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info == A_rows+1)
    py_error("Warning: matrix singular to working precision.");
  
  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }
  
  return X;
}



/////////////////////////////////////////////////////////////////////////////
//
// solve_svd(A,B)
//
/////////////////////////////////////////////////////////////////////////////

cMatrix solve_svd(const cMatrix& A, const cMatrix& B)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  // Solve system.

  return multiply(invert_svd(A), B);
}



/////////////////////////////////////////////////////////////////////////////
//
// solve_sym(A,B)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zsysv)
  (const char*,const int&,const int&,const Complex*,const int&,const int*,
   const Complex*,const int&,const Complex*,const int&,const int&);

cMatrix solve_sym(const cMatrix& A, const cMatrix& B)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  // Make copy of matrices, since LAPACK overwrites them.

  cMatrix A_LU(A_rows,A_cols,fortranArray);
  A_LU = A;
  
  cMatrix B_X(B_rows,B_cols,fortranArray);
  B_X = B;
  
  iVector P(A_rows,fortranArray); // pivot indices

  cVector work(2*A_rows,fortranArray);

  // Solve system.
  
  int info;
  
  F77NAME(zsysv)("U",A_rows,B_cols,A_LU.data(),A_rows,P.data(),B_X.data(),
                 B_rows,work.data(),2*A_rows,info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }

  return B_X;
}



/////////////////////////////////////////////////////////////////////////////
//
// solve_sym_x(A,B)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zsysvx)
  (const char*,const char*,const int&,const int&,const Complex*,const int&,
   const Complex*,const int&,const int*,const Complex*,const int&,
   const Complex*,const int&,Real&,const Real*,const Real*,const Complex*,
   const int&, const Real*,int&);

cMatrix solve_sym_x(const cMatrix& A, const cMatrix& B)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  // Make copy of matrices, since LAPACK overwrites them.

  cMatrix A_bis(A_rows,A_cols,fortranArray);
  A_bis = A;
  
  cMatrix B_bis(B_rows,B_cols,fortranArray);
  B_bis = B;

  // Prepare output arguments.
  
  cMatrix AF(A_rows,A_cols,fortranArray);  // factored A
  iVector P(A_rows,fortranArray);          // pivot indices
  cMatrix X(B_rows,B_cols,fortranArray);   // solution
  Real cond;                               // condition number
  rVector ferr(B_rows,fortranArray);       // forward error bound
  rVector berr(B_rows,fortranArray);       // backward error bound

  // Create workspace.
  
  cVector work (2*A_rows,fortranArray);
  rVector work2(2*A_rows,fortranArray);

  // Solve system.
  
  int info;

  F77NAME(zsysvx)("N","U",A_rows,B_cols,A_bis.data(),A_rows,AF.data(),A_rows,
                  P.data(),B_bis.data(),B_rows,X.data(),B_rows,cond,
                  ferr.data(),berr.data(),work.data(),2*A_rows,work2.data(),
                  info);

  // Diagnostic output:
  
  //cout << "Cond: " << cond << endl;
  //cout << "Ferr: " << ferr << endl;
  //cout << "Berr: " << berr << endl;

  if (info == 0)
    return X;
  
  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info == A_rows+1)
    py_error("Warning: matrix singular to working precision.");
  
  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }

  return X;
}



/////////////////////////////////////////////////////////////////////////////
//
// Computes eigenvalues and/or eigenvectors of matrix A.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgeev)
  (const char*,const char*,const int&,const Complex*,const int&,
   const Complex*,const Complex*,const int&,const Complex*,const int&,
   const Complex*,const int&,const Real*,int&);

cVector eigenvalues(const cMatrix& A, cMatrix* eigenvectors)
{
  // Check dimensions.

  const int N  = A.rows();

  if (N != A.columns())
  {
    py_error("Error: A matrix is not square.");
    exit (-1);
  }

  if (eigenvectors)
    eigenvectors->resize(N,N);
    
  // Make copy of matrix, since LAPACK overwrites them.

  cMatrix A_bis(N,N,fortranArray);
  A_bis = A;

  // Create workspace.

  const int work_size = 20; // variable to improve performance

  cVector eigenvalues(N,           fortranArray);
  cVector work       (work_size*N, fortranArray);
  rVector work2      (2*N,         fortranArray); // fixed

  // Calculate eigenvalues/vectors.
  
  int info;

  char op_code[] = " ";
  Complex* vectordata;
  
  if (eigenvectors)
  {
    op_code[0] = 'V';
    vectordata = eigenvectors->data();
  }
  else
  {
    op_code[0] = 'N';
    vectordata = NULL;
  }
  
  F77NAME(zgeev)("N",op_code,N,A_bis.data(),N,eigenvalues.data(),NULL,N,
                 vectordata,N,work.data(),work_size*N,work2.data(),info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: " << info << " eigenvalue(s) did not converge.";;
    py_error(s.str());
  }

  //cout << "Optimal work size: " << real(work(1))/N << "*N " << endl;
  
  return eigenvalues;
}



/////////////////////////////////////////////////////////////////////////////
//
// Computes eigenvalues and/or eigenvectors of matrix A. (expert driver)
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgeevx)
  (const char*,const char*,const char*,const char*,const int&,
   const Complex*,const int&,const Complex*,const Complex*,const int&,
   const Complex*,const int&,int&,int&,const Real*,Real&,const Real*,
   const Real*,const Complex*,const int&,const Real*,int&);

cVector eigenvalues_x(const cMatrix& A, cMatrix* eigenvectors)
{
  // Check dimensions.

  const int N  = A.rows();

  if (N != A.columns())
  {
    py_error("Error: A matrix is not square.");
    exit (-1);
  }

  if (eigenvectors)
    eigenvectors->resize(N,N);
    
  // Make copy of matrix, since LAPACK overwrites them.

  cMatrix A_bis(N,N,fortranArray);
  A_bis = A;
  
  // Prepare output arguments.

  int ilo, ihi;
  rVector scale(N,fortranArray);
  Real norm;
  rVector rconde(N,fortranArray);
  rVector rcondv(N,fortranArray);
  
  // Create workspace.

  const int work_size = 20; // variable to improve performance

  cVector eigenvalues(N,           fortranArray);
  cVector work       (work_size*N, fortranArray);
  rVector work2      (2*N,         fortranArray); // fixed

  // Calculate eigenvalues/vectors.
  
  int info;

  char op_code[] = " ";
  Complex* vectordata;
  
  if (eigenvectors)
  {
    op_code[0] = 'V';
    vectordata = eigenvectors->data();
  }
  else
  {
    op_code[0] = 'N';
    vectordata = NULL;
  }
  
  F77NAME(zgeevx)("B","N",op_code,"N",N,A_bis.data(),N,eigenvalues.data(),
                  NULL,N,vectordata,N,ilo,ihi,scale.data(),norm,rconde.data(),
                  rcondv.data(),work.data(),work_size*N,work2.data(),info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: " << info << " eigenvalue(s) did not converge.";
    py_error(s.str());
  }

  //cout << "Optimal work size: " << real(work(1))/N << "*N " << endl;
  
  return eigenvalues;
}



///////////////////////////////////////////////////////////////////////////// 
// 
// Computes eigenvalues and/or eigenvectors of the generalized 
// eigenproblem Ax = lambda Bx, where lambda = alpha / beta.
// 
///////////////////////////////////////////////////////////////////////////// 
 
extern "C" void F77NAME(zggev)
  (const char*,const char*,const int&,const Complex*,const int&,
   const Complex*,const int&,const Complex*,const Complex*,const Complex*,
   const int&,const Complex*,const int&,const Complex*,const int&,
   const Real*,int&);

void gen_eigenvalues(const cMatrix& A, const cMatrix& B,
                     cVector* alpha, cVector* beta,
                     cMatrix* eigenvectors) 
{
  // Check dimensions.

  int N  = A.rows();

  if (N != A.columns()) 
  {
    py_error("Error: A matrix is not square.");
    exit (-1); 
  }

  if (N != B.columns()) 
  {
    py_error("Error: A matrix has a different dimension than B.");
    exit (-1);
  }

  if (N != B.rows())
  {
    py_error("Error: B matrix is not square.");
    exit (-1);
  } 

  // Make copy of matrix, since LAPACK overwrites them.

  cMatrix A_bis(N,N,fortranArray); A_bis = A;
  cMatrix B_bis(N,N,fortranArray); B_bis = B;

  // Make sure input vectors and matrices are the correct size.

  alpha->resize(N);
   beta->resize(N);

  if (eigenvectors)
    eigenvectors->resize(N,N);

  // Create workspace.

  const int work_size = 33; // Variable to improve performance. 
  cVector work (work_size*N,fortranArray); 
  rVector work2(8*N,fortranArray); // Fixed. 

  // Calculate eigenvalues/vectors. 

  int info; 
  char op_code[] = " "; 
  Complex* vectordata; 

  if (eigenvectors) 
  {
    op_code[0] = 'V';
    vectordata = eigenvectors->data();
  } 
  else 
  {
    op_code[0] = 'N';
    vectordata = NULL; 
  }

  F77NAME(zggev)("N",op_code,N,A_bis.data(),N,B_bis.data(),N,
                 alpha->data(),beta->data(),NULL,N,vectordata,
                 N,work.data(),work_size*N,work2.data(),info); 

  if (info < 0) 
  {
    std::ostringstream s; 
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1); 
  } 

  if ((info > 0) && (info <= N))
  {
    std::ostringstream s;
    s << "The QZ iteration failed. " 
      << "No eigenvectors have been calculated," << std::endl
      << "but alpha(j) and beta(j) should be "
      << "correct for j=" << info+1 << ",..., " << N; 
    py_error(s.str());
  }

  if (info > N) 
  {
    std::ostringstream s;
    s << "Warning: " << info << " eigenvalue(s) did not converge.";
    py_error(s.str());
  }

  //cout << "Optimal work size: " << real(work(1))/N << "*N " << endl;
}



/////////////////////////////////////////////////////////////////////////////
//
// Computes the SVD decomposition U.Sigma.Vh of matrix A.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgesvd)
  (const char*,const char*,const int&,const int&,const Complex*,const int&,
   const Real*,const Complex*,const int&,const Complex*,const int&,
   const Complex*,const int&,const Real*,int&);

rVector svd(const cMatrix& A, cMatrix* Vh, cMatrix* U)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();

  if (Vh)
    Vh->resize(A_cols,A_cols);
  if (U)
    U->resize(A_rows,A_rows);
    
  // Make copy of matrix, since LAPACK overwrites them.

  cMatrix A_bis(A_rows,A_cols,fortranArray);
  A_bis = A;

  // Create workspace.

  const int min_dim = (A_rows < A_cols) ? A_rows : A_cols;
  const int max_dim = (A_rows > A_cols) ? A_rows : A_cols;
  
  const int work_size = 20; // variable to improve performance
  
  rVector SVD  (min_dim,           fortranArray);
  cVector work (work_size*max_dim, fortranArray); 
  rVector work2(5*max_dim,         fortranArray); // fixed

  // Calculate SVD decomposition.
  
  int info;

  char Vh_code[] = " ";
  Complex* Vh_data;
  
  if (Vh)
  {
    Vh_code[0] = 'A';
    Vh_data = Vh->data();
  }
  else
  {
    Vh_code[0] = 'N';
    Vh_data = NULL;
  }

  char U_code[] = " ";
  Complex* U_data;
  
  if (U)
  {
    U_code[0] = 'A';
    U_data = U->data();
  }
  else
  {
    U_code[0] = 'N';
    U_data = NULL;
  }
  
  F77NAME(zgesvd)(U_code,Vh_code,A_rows,A_cols,A_bis.data(),A_rows,SVD.data(),
                  U_data,A_rows,Vh_data,A_cols,work.data(),work_size*max_dim,
                  work2.data(),info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: element " << info << " did not converge.";
    py_error(s.str());
  }

  //cout << "Optimal work size: " << real(work(1))/max_dim << "*N " << endl;
  
  return SVD;
}


  
/////////////////////////////////////////////////////////////////////////////
//
// invert(A)
//  
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgetri)
  (const int&,const Complex*,const int&,const int*,const Complex*,
   const int&,int&);

cMatrix invert(const cMatrix& A)
{
  // Check dimensions.

  const int N = A.rows();
  
  if (N != A.columns())
  {
    py_error("Can only invert square matrices.");
    exit (-1);
  }

  // Create workspace.

  const int work_size = 20; // variable to improve performance
  
  cMatrix LU_invA(N,N,         fortranArray);
  iVector P      (N,           fortranArray);
  cVector work   (work_size*N, fortranArray);
  
  // Calculate inverse.

  int info;

  LU(A,&LU_invA,&P);

  F77NAME(zgetri)(N,LU_invA.data(),N,P.data(),work.data(),work_size*N,info);
  
  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }

  //cout << "Optimal work size: " << real(work(1))/N << "*N " << endl;
  
  return LU_invA;
}



/////////////////////////////////////////////////////////////////////////////
//
// invert_x(A)
//  
/////////////////////////////////////////////////////////////////////////////

cMatrix invert_x(const cMatrix& A)
{
  // Check dimensions.

  const int N = A.rows();
  
  if (N != A.columns())
  {
    py_error("Can only invert square matrices.");
    exit (-1);
  }

  // Create unity matrix.

  cMatrix U1(N,N,fortranArray);
  U1 = 0.0;
  for (int i=1; i<=N; i++)
    U1(i,i) = 1.0;

  // solve_x on system

  return solve_x(A, U1);
}



/////////////////////////////////////////////////////////////////////////////
//
// invert_svd(A)
//  
/////////////////////////////////////////////////////////////////////////////

cMatrix invert_svd(const cMatrix& A)
{ 
  // Check dimensions.

  const int N = A.rows();
  
  if (N != A.columns())
  {
    py_error("Can only invert square matrices.");
    exit (-1);
  }

  // Do SVD of A: A=U.sigma.Vh

  rVector sigma(N,    fortranArray);
  cMatrix     U(N, N, fortranArray);
  cMatrix    Vh(N, N, fortranArray);
  
  sigma = svd(A, &Vh, &U);
  
  // Invert sigma, but set 1/0 to 0

  int singularities = 0;
  
  for (int i=1; i<=N; i++)
    if ( abs(sigma(i)) <= 1e-11 )
    {
      singularities++;
      sigma(i) = 0;
    }
    else
      sigma(i) = 1./sigma(i);

  if (singularities)
  {
    std::ostringstream s;
    s << singularities << "/" << N << " zero singular values encountered.";
    py_error(s.str());
  }

  // hermitian conjugate U and Vh

  herm_conj_self(&U);
  herm_conj_self(&Vh);

  // inverse(A) = V.1/sigma.Uh
  
  blitz::firstIndex i;
  blitz::secondIndex j;
  
  U = sigma(i)*U(i,j);

  return multiply(Vh,U);
}



/////////////////////////////////////////////////////////////////////////////
//
// determinant
//
/////////////////////////////////////////////////////////////////////////////

Complex determinant(const cMatrix& A)
{
  // Check dimensions.

  const int N = A.rows();
  
  if (N != A.columns())
  {
    py_error("Only square matrices currently supported.");
    exit (-1);
  }

  // Do LU decomposition.

  cMatrix LU_(N,N,fortranArray);
  iVector P (N,  fortranArray);
  LU(A, &LU_, &P);

  // Calculate determinant.

  Complex determinant = 1.0;
  for (int i=1; i<=N; i++)
    determinant *= LU_(i,i);

  // Calcute sign change due to permutations.

  int changed = 0;
  for (int i=1; i<=N; i++)
    if (P(i) != i)
      changed++;
  
  int swaps = int(changed/2);
  bool odd_swaps = 2*int(swaps/2) != swaps;
  
  return odd_swaps ? -determinant : determinant;
}



/////////////////////////////////////////////////////////////////////////////
//
// Do LU decomposition of matrix A.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgetrf)
  (const int&,const int&,const Complex*,const int&,const int*,int&);

void LU(const cMatrix& A, cMatrix* LU, iVector* P)
{
  // Check dimensions.

  const int A_rows = A.rows();
  const int A_cols = A.columns();

  LU->resize(A_rows,A_cols);
  P->resize( (A_rows < A_cols) ? A_rows : A_cols );
    
  // Make copy of matrix, since LAPACK overwrites them.
  
  *LU = A;

  // Calculate LU decomposition.
  
  int info;

  F77NAME(zgetrf)(A_rows,A_cols,LU->data(),A_rows,P->data(),info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }

  if (info > 0)
  {
    std::ostringstream s;
    s << "Warning: singular system: U(" << info << "," << info << ") is zero.";
    py_error(s.str());
  }
}




/////////////////////////////////////////////////////////////////////////////
//
// Solve linear system op(A)*X=B starting from LU decomposition of A.
//
/////////////////////////////////////////////////////////////////////////////

extern "C" void F77NAME(zgetrs)
  (const char*,const int&,const int&,const Complex*,const int&,
   const int*,const Complex*,const int&,int&);

cMatrix LU_solve(const cMatrix& LU, const iVector& P,
                 const cMatrix& B, Op op)
{ 
  // Check dimensions.

  const int A_rows = LU.rows();
  const int A_cols = LU.columns();
  const int B_rows = B.rows();
  const int B_cols = B.columns();

  if (A_rows != A_cols)
  {
    py_error("Error: system matrix is not square.");
    exit (-1);
  }

  if (A_rows != B_rows)
  {
    py_error("Error: dimension of rhs matrix does not match.");
    exit (-1);
  }

  if (P.rows() != A_rows)
  {
    py_error("Error: incorrect dimension of permution matrix.");
    exit (-1);
  }

  // Make copy of matrices, since LAPACK overwrites them.
  
  cMatrix B_X(B_rows,B_cols,fortranArray);
  B_X = B;

  // Solve system.
  
  int info;

  char op_code[] = "N";
  if (op==transp)
    op_code[0] = 'T';
  if (op==herm)
    op_code[0] = 'C';
  
  F77NAME(zgetrs)(op_code,A_rows,B_cols,LU.data(),A_rows,P.data(),
                  B_X.data(),B_rows,info);

  if (info < 0)
  {
    std::ostringstream s;
    s << "Error: bad value for argument " << -info;
    py_error(s.str());
    exit (-1);
  }
  
  return B_X;
}

  






