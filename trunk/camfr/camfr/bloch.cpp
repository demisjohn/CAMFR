
/////////////////////////////////////////////////////////////////////////////
//
// File:     bloch.cpp
// Authors:  Peter.Bienstman@rug.ac.be, Lieven.Vanholme@rug.ac.be
// Date:     20020207
// Version:  2.0
//
// Copyright (C) 2000-2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "bloch.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// BlochStack::BlochStack
//
//   The factor 1 in 1*e ensures that incidence and exit media of the
//   period are the same.
//
/////////////////////////////////////////////////////////////////////////////

BlochStack::BlochStack(const Expression& e) : stack(Expression(1*e)) {}



/////////////////////////////////////////////////////////////////////////////
//
// BlochStack::find_modes
//
/////////////////////////////////////////////////////////////////////////////

void BlochStack::find_modes()
{
  if (stack.get_expression().all_layers_uniform() && global.orthogonal)
    find_modes_diag();
  else if (global.bloch_calc == GEV)
    find_modes_GEV();
  else
    find_modes_T();
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochStack::find_modes_T
//
/////////////////////////////////////////////////////////////////////////////

void BlochStack::find_modes_T()
{  
  stack.calcRT();

  // Calculate submatrices.
  
  const int N = global.N;

  cMatrix A(N,N,fortranArray); cMatrix B(N,N,fortranArray);
  cMatrix C(N,N,fortranArray); cMatrix D(N,N,fortranArray);

  const cMatrix& R12(stack.as_multi()->get_R12());
  const cMatrix& R21(stack.as_multi()->get_R21());
  const cMatrix& T12(stack.as_multi()->get_T12());
  const cMatrix& T21(stack.as_multi()->get_T21());

  cMatrix inv_T21    (N,N,fortranArray);
  cMatrix inv_T21_R12(N,N,fortranArray);

  if (global.stability != SVD)
    inv_T21.reference(invert    (T21));
  else
    inv_T21.reference(invert_svd(T21));
  
  inv_T21_R12.reference(multiply(inv_T21, R12));

  A =         T12 - multiply(R21, inv_T21_R12);
  B.reference(      multiply(R21, inv_T21));
  C =                        -inv_T21_R12;
  D.reference(                    inv_T21);
  
  // Create matrix for eigenvalue problem.

  cMatrix E(2*N,2*N,fortranArray);
  
  blitz::Range r1(1,N); blitz::Range r2(N+1,2*N);

  E(r1,r1) = A; E(r1,r2) = B;
  E(r2,r1) = C; E(r2,r2) = D;

  // Solve eigenproblem.

  cVector alpha(2*N,fortranArray);
  cMatrix FB(2*N,2*N,fortranArray);

  alpha = eigenvalues(E, &FB);
  
  // Create modeset.

  modeset.clear();
  for (int mode=1; mode<=2*N; mode++)
  {
    Complex beta = log(alpha(mode))/(-I*stack.get_total_thickness());

    if (abs(alpha(mode)) < 1e-10)
        beta = 0;
    
    cVector F(N,fortranArray), B(N,fortranArray);
    for (int i=1; i<=N; i++)
    {
      F(i) = FB(  i,mode);
      B(i) = FB(N+i,mode);
    }
    
    Polarisation pol = stack.get_inc()->get_mode(1)->pol;

    BlochMode* blochmode = new BlochMode(pol,beta,&stack,F,B);

    modeset.push_back(blochmode);
  }

  sort_modes_bloch();
}



/////////////////////////////////////////////////////////////////////////////
// 
// BlochStack::find_modes_GEV 
// 
//  Finds Blochmodes by solving generalised eigenproblem:
//
//  | T21  R21 | |  F | = e |  O   I  | |  F | 
//  |  O    I  | |e B |     | R12 T21 | |e B | 
// 
///////////////////////////////////////////////////////////////////////////// 

void BlochStack::find_modes_GEV() 
{ 
  // Set up matrices.

  stack.calcRT();

  const cMatrix& R12(stack.as_multi()->get_R12()); 
  const cMatrix& R21(stack.as_multi()->get_R21()); 
  const cMatrix& T12(stack.as_multi()->get_T12()); 
  const cMatrix& T21(stack.as_multi()->get_T21()); 

  const int N = global.N;
  
  cMatrix U1(N,N,fortranArray); 
  U1 = 0.0;
  for(int i=1; i<=N; i++) 
    U1(i,i) = 1.0;

  // Solve generalised eigenvalue problem.

  cMatrix A(2*N,2*N,fortranArray); 
  cMatrix B(2*N,2*N,fortranArray); 

  blitz::Range r1(1,N); blitz::Range r2(N+1,2*N); 

  A(r1,r1) = T12; A(r1,r2) = R21; 
  A(r2,r1) = 0.0; A(r2,r2) =  U1;

  B(r1,r1) =  U1; B(r1,r2) = 0.0; 
  B(r2,r1) = R12; B(r2,r2) = T21; 

  cVector a(2*N,fortranArray); // beta = a/b
  cVector b(2*N,fortranArray); 
  cMatrix Z(2*N,2*N,fortranArray); // Holds F and e B.

  gen_eigenvalues(A, B, &a, &b, &Z); 

  // Create modeset. 

  modeset.clear(); 

  for (int mode=1; mode<=2*N; mode++) 
  {
    // Determine propagation constant.
    
    Complex beta;

    if (abs(b(mode)) < 1e-30) // Set infinite eigenvalue to zero.
      beta = 0.0;
    else
      beta = I*log(a(mode)/b(mode)) / stack.get_total_thickness();

    // Determine fields.

    cVector F(N,fortranArray); F = Z(r1,mode);
    cVector B(N,fortranArray); B = Z(r2,mode) / (a(mode)/b(mode));

    if (abs(beta) < 1e-20) // Non-converged result.
    {
      F = 0.0; 
      B = 0.0;
    }
    
    // Create mode.

    Polarisation pol = stack.get_inc()->get_mode(1)->pol; 

    BlochMode* blochmode = new BlochMode(pol,beta,&stack,F,B); 

    modeset.push_back(blochmode);
  } 

  sort_modes_bloch(); 
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochStack::find_modes_diag
//
/////////////////////////////////////////////////////////////////////////////

void BlochStack::find_modes_diag()
{
  stack.calcRT();
  
  modeset.clear();

  for (int i=1; i<=global.N; i++)
  {
    // Create matrix for eigenvalue problem.

    DiagStack* s = dynamic_cast<DiagStack*>(stack.get_sc());

    if (!s)
    {
      py_error("Error: invalid diagonal stack.");
      exit (-1);
    }
    
    const Complex& R12(s->get_diag_R12()(i));
    const Complex& R21(s->get_diag_R21()(i));
    const Complex& T12(s->get_diag_T12()(i));
    const Complex& T21(s->get_diag_T21()(i));

    cMatrix E(2,2,fortranArray);

    E(1,1) = T12 - R21 * R12 / T21;
    E(1,2) =       R21       / T21;
    E(2,1) =           - R12 / T21;
    E(2,2) =             1.0 / T21;

    if (abs(T21) < 1e-10)
    {
      //cout << "Error: mode near cut-off. Impossible to calculate "
      //     << "bloch modes" << endl;
      E = 0;
    }
    
    // Solve eigenproblem.

    cVector alpha(2,fortranArray);
    cMatrix FB(2,2,fortranArray);

    alpha = eigenvalues(E, &FB);

    // Create modeset.

    for (int mode=1; mode<=2; mode++)
    {
      Complex beta = log(alpha(mode))/(-I*stack.get_total_thickness());

      if (abs(alpha(mode)) < 1e-10)
        beta = 0;

      cVector F(global.N,fortranArray); F = 0; F(i) = FB(1,mode);
      cVector B(global.N,fortranArray); B = 0; B(i) = FB(2,mode);

      Polarisation pol = stack.get_inc()->get_mode(1)->pol;

      BlochMode* blochmode = new BlochMode(pol,beta,&stack,F,B);

      modeset.push_back(blochmode);
    }
  }

  sort_modes_bloch();
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochStack::get_beta_vector
//
/////////////////////////////////////////////////////////////////////////////

cVector BlochStack::get_beta_vector() const
{
  cVector beta(N(), fortranArray);

  for (int i=1; i<=N(); i++)
    beta(i) = get_mode(i)->get_kz();

  return beta;
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochMode::BlochMode
//
/////////////////////////////////////////////////////////////////////////////

BlochMode::BlochMode(const Polarisation pol, const Complex& kz, Stack* s,
                     cVector& F, cVector& B)
  : Mode(pol,kz,-kz), geom(s)
{
  FieldExpansion f(*(geom->get_inc()), F, B);
  interface_field.push_back(f);
}



/////////////////////////////////////////////////////////////////////////////
//
// BlochMode::field
//
/////////////////////////////////////////////////////////////////////////////

Field BlochMode::field(const Coord& coord) const
{
  Coord coord2(coord);
  
  if (real(coord.z) > real(geom->get_total_thickness()))
  {
    py_print("Warning: z-value out of unit cell. Truncating z.");
    coord2.z = geom->get_total_thickness();
  }

  geom->set_interface_field(interface_field);

  if (global.field_calc == S_T)
    py_print("Warning: field_calc S_T can give unstable results.");

  if ( (global.field_calc == S_S) && (interface_field.size() <= 1) )
  {
    Complex d = geom->get_total_thickness();
    cVector inc_bw(global.N, fortranArray);
    inc_bw = interface_field[0].bw * exp(-I*kz*d);
 
    geom->set_inc_field(interface_field[0].fw, &inc_bw);
  }

  Field f(geom->field(coord2));

  if (interface_field.size() <= 1)
    geom->get_interface_field(&interface_field);

  return f;
}
