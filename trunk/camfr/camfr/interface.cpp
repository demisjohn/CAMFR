
/////////////////////////////////////////////////////////////////////////////
//
// File:     interface.cpp
// Author:   Peter.Bienstman@rug.ac.be
// Date:     19991105
// Version:  1.0
//
// Copyright (C) 1999 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "interface.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::get_materials
//
/////////////////////////////////////////////////////////////////////////////

vector<Material*> DenseInterface::get_materials() const
{
  vector<Material*> inc_materials = inc->get_materials();
  vector<Material*> ext_materials = ext->get_materials();

  for (unsigned int i=0; i<ext_materials.size(); i++)
    inc_materials.push_back(ext_materials[i]);
  
  return inc_materials;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::calcRT
//
/////////////////////////////////////////////////////////////////////////////

void DenseInterface::calcRT()
{  
  if (!recalc_needed())
    return;

  allocRT();

  inc->find_modes();
  ext->find_modes();

  if (global.orthogonal == false)
    (global.stability == normal) ? calcRT_non_orth_fast()
                                 : calcRT_non_orth_safe();
  else
    (global.stability == normal) ? calcRT_fast()
                                 : calcRT_safe();

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
  last_slab_ky = global.slab_ky;
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::calcRT_safe
//
//  Uses row/column equilibration when solving the system. Rarely needed.
//
/////////////////////////////////////////////////////////////////////////////

void DenseInterface::calcRT_safe()
{
  // Set constants and calculate overlap matrices.

  blitz::firstIndex  i;
  blitz::secondIndex j;

  const int N = global.N;

  cMatrix O_I_II(N,N,fortranArray);
  cMatrix O_II_I(N,N,fortranArray);

  cMatrix A(N,N,fortranArray);
  cMatrix B(N,N,fortranArray);
  cMatrix C(N,N,fortranArray);

  MultiWaveguide* w1 = dynamic_cast<MultiWaveguide*>(inc);
  MultiWaveguide* w2 = dynamic_cast<MultiWaveguide*>(ext);  

  w1->calc_overlap_matrices(w2, &O_I_II, &O_II_I);
  
  //
  // Calculate R12 and T12.
  //

  // Calculate system matrix A from overlap matrices.
  
  A = O_I_II(i,j) + O_II_I(j,i);
  
  // Calculate diagonal rhs matrix B.
  
  for (int i=1; i<=N; i++)
    for (int p=1; p<=N; p++)
      B(i,p) = (i!=p) ? 0.0 : 2.0; // 2.0 * O_I_I(p);
  
  // Solve system for T12.

  if (global.stability == extra)
    T12.reference(  solve_x(A,B));
  else
    T12.reference(solve_svd(A,B));    
  
  // Calculate R12 from T12.

  C = O_II_I(j,i) - O_I_II(i,j);
  
  R12.reference(multiply(C,T12));

  R12 /= Complex(2.0); // 2.0 * O_I_I(i);
  
  //
  // Calculate R21 and T21.
  //

  // Calculate system matrix A from overlap matrices.
  
  A = O_II_I(i,j) + O_I_II(j,i);
  
  // Calculate diagonal rhs matrix B.
  
  for (int i=1; i<=N; i++)
    for (int p=1; p<=N; p++)
      B(i,p) = (i!=p) ? 0.0 : 2.0; // 2.0 * O_II_II(p);
  
  // Solve system for T21.

  if (global.stability == extra)
    T21.reference(  solve_x(A,B));
  else
    T21.reference(solve_svd(A,B));

  // Calculate R21 from T21.

  C = O_I_II(j,i) - O_II_I(i,j);
  
  R21.reference(multiply(C,T21));

  R21 /= Complex(2.0); // 2.0 * O_II_II(i); 
}


  
/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::calcRT_fast
//
//   Doesn't do row/column equilibration, which allows for some shortcuts.
//
/////////////////////////////////////////////////////////////////////////////

void DenseInterface::calcRT_fast()
{ 
  // Set constants and calculate overlap matrices.

  blitz::firstIndex  i;
  blitz::secondIndex j;

  const int N = global.N;

  cMatrix O_I_II(N,N,fortranArray);
  cMatrix O_II_I(N,N,fortranArray);

  cMatrix A(N,N,fortranArray);
  cMatrix B(N,N,fortranArray);
  cMatrix C(N,N,fortranArray);

  MultiWaveguide* w1 = dynamic_cast<MultiWaveguide*>(inc);
  MultiWaveguide* w2 = dynamic_cast<MultiWaveguide*>(ext);  

  w1->calc_overlap_matrices(w2, &O_I_II, &O_II_I);
  
  //
  // Calculate R12 and T12.
  //

  // Calculate system matrix A from overlap matrices.
  
  A = O_I_II(i,j) + O_II_I(j,i);
  
  // Calculate diagonal rhs matrix B.
  
  for (int i=1; i<=N; i++)
    for (int p=1; p<=N; p++)
      B(i,p) = (i!=p) ? 0.0 : 2.0; // 2.0 * O_I_I(p);
  
  // Solve system for T12, but save LU decomposition for later use.

  cMatrix lu(N,N,fortranArray);
  iVector  p(  N,fortranArray);
  
  LU(A,&lu,&p);
  
  T12.reference(LU_solve(lu,p,B));
  
  // Calculate R12 from T12.

  C = O_II_I(j,i) - O_I_II(i,j);  
  
  R12.reference(multiply(C,T12));

  R12 /= Complex(2.0); // 2.0 * O_I_I(i);
  
  //
  // Calculate R21 and T21, but use some shortcuts.
  //

  // The new A matrix is the transpose of the old A matrix.
  
  // Calculate diagonal rhs matrix B.
  // The off-diagonal elements were already set to zero earlier.

  for (int i=1; i<=N; i++)
    B(i,i) = 2.0; // 2.0 * O_II_II(i);
  
  // Solve system for T21, but reuse LU decomposition of A.
  // Note: the transposition of A enters here.

  T21.reference(LU_solve(lu,p,B,transp));

  // Calculate R21 from T21.
  // Note: the new C is minus the transpose of the old C.
  
  R21.reference(multiply(C,T21,transp)); // (transpose)

  R21 /= Complex(-2.0); // -2.0 * 0_II_II(i) (minus);
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::calcRT_non_orth_safe
//
//  Calculates RT when the modes are not orthogonal.
//  Uses row/column equilibration when solving the system. Rarely needed.
//
/////////////////////////////////////////////////////////////////////////////

void DenseInterface::calcRT_non_orth_safe()
{
  // Set constants and calculate overlap matrices.

  const int N = global.N;

  blitz::Range r1(1,N); blitz::Range r2(N+1,2*N);

  cMatrix O_I_II (N,N,fortranArray);
  cMatrix O_II_I (N,N,fortranArray);
  cMatrix O_I_I  (N,N,fortranArray);
  cMatrix O_II_II(N,N,fortranArray);

  cMatrix A(2*N,2*N,fortranArray);
  cMatrix B(2*N,  N,fortranArray);
  cMatrix X(2*N,  N,fortranArray);

  MultiWaveguide* w1 = dynamic_cast<MultiWaveguide*>(inc);
  MultiWaveguide* w2 = dynamic_cast<MultiWaveguide*>(ext);  

  w1->calc_overlap_matrices(w2, &O_I_II, &O_II_I, &O_I_I, &O_II_II);
  
  //
  // Calculate R12 and T12.
  //

  // Calculate system matrix A from overlap matrices.

  A(r1,r1) = transpose(O_II_I); A(r1,r2) = -transpose(O_I_I);
  A(r2,r1) =           O_I_II;  A(r2,r2) =            O_I_I;
  
  // Calculate rhs matrix B.

  B(r1,r1) = transpose(O_I_I);
  B(r2,r1) =           O_I_I;
  
  // Solve system for R12 and T12.

  if (global.stability == extra)
    X.reference(solve_x(A,B));
  else
    X.reference(solve_svd(A,B));

  T12 = X(r1,r1);
  R12 = X(r2,r1);
  
  //
  // Calculate R21 and T21.
  //

  // Calculate system matrix A from overlap matrices.

  A(r1,r1) = transpose(O_I_II); A(r1,r2) = -transpose(O_II_II);
  A(r2,r1) =           O_II_I;  A(r2,r2) =            O_II_II;
  
  // Calculate rhs matrix B.

  B(r1,r1) = transpose(O_II_II);
  B(r2,r1) =           O_II_II;
  
  // Solve system for R12 and T12
  
  if (global.stability == extra)
    X.reference(solve_x(A,B));
  else
    X.reference(solve_svd(A,B));

  T21 = X(r1,r1);
  R21 = X(r2,r1);
}



/////////////////////////////////////////////////////////////////////////////
//
// DenseInterface::calcRT_non_orth_fast
//
//  Calculates RT when the modes are not orthogonal.
//  Doesn't do row/column equilibration and uses alternative formulation,
//  which allows for some shortcuts exploiting the symmetry.
//
/////////////////////////////////////////////////////////////////////////////

void DenseInterface::calcRT_non_orth_fast()
{
  // Set constants and calculate overlap matrices.

  const int N = global.N;

  blitz::Range r1(1,N); blitz::Range r2(N+1,2*N);

  cMatrix O_I_II (N,N,fortranArray);
  cMatrix O_II_I (N,N,fortranArray);
  cMatrix O_I_I  (N,N,fortranArray);
  cMatrix O_II_II(N,N,fortranArray);

  cMatrix A(2*N,2*N,fortranArray);
  cMatrix B(2*N,2*N,fortranArray);
  cMatrix X(2*N,2*N,fortranArray);

  MultiWaveguide* w1 = dynamic_cast<MultiWaveguide*>(inc);
  MultiWaveguide* w2 = dynamic_cast<MultiWaveguide*>(ext);  

  w1->calc_overlap_matrices(w2, &O_I_II, &O_II_I, &O_I_I, &O_II_II);

  // Calculate system matrix A.

  A(r1,r1) = transpose(O_II_II); A(r1,r2) = -transpose(O_I_II);
  A(r2,r1) =           O_I_II;   A(r2,r2) =            O_I_I;
  
  // Calculate rhs matrix B.

  B(r1,r1) = transpose(O_I_II);  B(r1,r2) = -transpose(O_II_II);
  B(r2,r1) =           O_I_I;    B(r2,r2) =            O_I_II;
  
  // Solve system.

  X.reference(solve(A,B));

  T12 = X(r1,r1); R21 = X(r1,r2);
  R12 = X(r2,r1); T21 = X(r2,r2);
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagInterface::get_materials
//
/////////////////////////////////////////////////////////////////////////////

vector<Material*> DiagInterface::get_materials() const
{
  vector<Material*> inc_materials = inc->get_materials();
  vector<Material*> ext_materials = ext->get_materials();

  for (unsigned int i=0; i<ext_materials.size(); i++)
    inc_materials.push_back(ext_materials[i]);
  
  return inc_materials;
}



/////////////////////////////////////////////////////////////////////////////
//
// calc_RT_fresnel
//
//   Calculate R and T for matched mode m1 in medium inc and mode m2 in
//   medium ext.
//
//   Note: the conventions obeyed in these calculations are such that the
//   formulas can be used both for planar stratified and circular media,
//   and that they give the same results as the overlap integral
//   calculation method. More specifically:
//
//   In order to get the same results as with the general overlap integral
//   formulas, the convention for the sign of r_tm is such that r_te = r_tm
//   for theta=0.
//
//   The 'T' factors in transmission are needed to make the formulas
//   applicable for different geometries (planar infinite, homogeneously
//   filled metal cilinder, ..). The values of E_cst and H_cst are supposed
//   to be set correctly by the modes for this to work.
//  
/////////////////////////////////////////////////////////////////////////////

void calc_RT_fresnel(const Mode* m1,       const Mode* m2,
                     const Waveguide* inc, const Waveguide* ext,
                     Complex* R12,         Complex* R21,
                     Complex* T12,         Complex* T21)
{
  Complex a, T;
  Real sign;

  //
  // Normal formulas, accurate for non-grazing incidence.
  //

  const Complex kz1 = m1->get_kz();
  const Complex kz2 = m2->get_kz();

  if ( abs(kz1) > abs(kz2) )
  {
    if ( (abs(kz1) < 1e-10) && (abs(kz2) < 1e-10) )
      a = 1.0;
    else
      a = kz2 / kz1;
    
    if (m1->pol == TE)
    {
      a *=    inc->get_core()->mu()
           /  ext->get_core()->mu();

      T =   ( inc->get_core()->mu() * m1->H_cst() )
          / ( ext->get_core()->mu() * m2->H_cst() );
      
      sign = 1.0;
    }
    else // TM
    {
      a *=    inc->get_core()->eps()
           /  ext->get_core()->eps();

      T =   ( inc->get_core()->eps() * m1->E_cst() )
          / ( ext->get_core()->eps() * m2->E_cst() );

      sign = -1.0; // Make r_te = t_tm for theta=0.0.
    }

    *R12 = sign * (1.0 -  a ) / (1.0 +  a );
    *R21 = sign * ( a  - 1.0) / ( a  + 1.0);
    *T12 =            2.0     / (1.0 +  a ) * T;
    *T21 =          2.0 * a   / ( a  + 1.0) / T;
  }

  //
  // Formulas accurate for grazing incidence.
  //

  else
  {
    if ( (abs(kz1) < 1e-10) && (abs(kz2) < 1e-10) )
      a = 1.0;
    else
      a = kz1 / kz2;
    
    if (m1->pol == TE)
    {
      a *=    ext->get_core()->mu()
           /  inc->get_core()->mu();

      T =   ( inc->get_core()->mu() * m1->H_cst() )
          / ( ext->get_core()->mu() * m2->H_cst() );

      sign = 1.0;
    }
    else // TM
    {
      a *=    ext->get_core()->eps()
           /  inc->get_core()->eps();

      T =   ( inc->get_core()->eps() * m1->E_cst() )
          / ( ext->get_core()->eps() * m2->E_cst() );
      
      sign = -1.0; // Make r_te = t_tm for theta=0.0.
    }
    
    *R12 = sign * ( a  - 1.0) / ( a  + 1.0);
    *R21 = sign * (1.0 -  a ) / (1.0 +  a );
    *T12 =          2.0 * a   / ( a  + 1.0) * T;
    *T21 =            2.0     / (1.0 +  a ) / T;
  }
}



/////////////////////////////////////////////////////////////////////////////
//
// DiagInterface::calcRT
//  
/////////////////////////////////////////////////////////////////////////////

void DiagInterface::calcRT()
{
  if (!recalc_needed())
    return;

  allocRT();

  if (abs(inc->c1_size() - ext->c1_size()) > 1e-6)
  {
    std::ostringstream s;
    s << "Warning: complex widths don't match: "
      << inc->c1_size() << " and " << ext->c1_size();
    py_error(s.str());
    exit (-1);
  }

  inc->find_modes();
  ext->find_modes();

  if (!global.orthogonal)
  {
    //py_error("Error: uniform waveguide with PC_Walls is not diagonal.");
    //exit (-1);
  }
  
  // Calc RT according to Fresnel laws.
  
  for (int i=1; i<=int(global.N); i++)
    calc_RT_fresnel(inc->get_mode(i), ext->get_mode(i), inc, ext,
                    &R12(i), &R21(i), &T12(i), &T21(i));

  // Remember wavelength and gain these matrices were calculated for.

  last_lambda = global.lambda;
  if (global.gain_mat)
    last_gain_mat = *global.gain_mat;
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoInterface::get_materials
//
/////////////////////////////////////////////////////////////////////////////

vector<Material*> MonoInterface::get_materials() const
{
  vector<Material*> materials;
  materials.push_back(inc->get_core());
  materials.push_back(ext->get_core());
  return materials;
}



/////////////////////////////////////////////////////////////////////////////
//
// MonoInterface::calcRT
//
//   Note: the convention for the sign of r_tm is such that
//   r_te = r_tm for theta=0.
//  
/////////////////////////////////////////////////////////////////////////////

void MonoInterface::calcRT()
{ 
  inc->find_modes();
  ext->find_modes();

  calc_RT_fresnel(inc->get_mode(1), ext->get_mode(1), inc, ext,
                  &R12, &R21, &T12, &T21);
}
