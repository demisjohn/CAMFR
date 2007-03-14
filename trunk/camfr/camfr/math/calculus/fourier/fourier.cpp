/////////////////////////////////////////////////////////////////////////////
//
// File:     fourier.cpp
// Author:   Peter.Bienstman@UGent.be, Peter.Debackere@intec.UGent.be
// Date:     20050208
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include "fourier.h"

using std::vector;

/////////////////////////////////////////////////////////////////////////////
//
// fourier
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier(const vector<Complex>& f, const vector<Complex>& disc, 
                int M, const Complex* d, bool extend)
{
  Complex D = d ? *d : disc.back() - disc.front();
  if (extend)
    D *= 2.0;
  Complex K = 2.*pi/D;

  cVector result(2*M+1,fortranArray);

  for (int m=-M; m<=M; m++)
  {
    Complex result_m = 0.0;

    for (unsigned int k=0; k<int(disc.size()-1); k++)
    { 
      Complex factor;
      if (m==0)
      { 
        factor = disc[k+1]-disc[k];
        if (extend)
          factor *= 2.0;
      }
      else
      {
        const Complex t = -I*Real(m)*K;
        factor = (exp(t*disc[k+1])-exp(t*disc[k])) / t;
        if (extend)
          factor += (exp(-t*disc[k])-exp(-t*disc[k+1])) / t;
      }

      result_m += factor * f[k];
    }
    
    result(m+M+1) = result_m / D;
  }

  if (abs(D) < 1e-12)
    result = 0.0;
  
  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_ASR
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier_ASR(const vector<Complex>& f,
                    const vector<Complex>& disc, 
                    const vector<Complex>& udisc,
                    int M, bool extend, const Real eta)
{ 
  if (0)
  {
    std::cout << "Fourier_ASR" << std::endl;
    std::cout << "Original discontinuity profile" << std::endl;
    for (unsigned int k=0; k<disc.size(); k++)
      std::cout << "discontinuity " << k+1 << " " << disc[k] << std::endl;

    std::cout << "New discontinuity profile" << std::endl;
    for (unsigned int k=0; k<udisc.size(); k++)
      std::cout << "discontinuity " << k+1 << " " << udisc[k] <<std::endl;
  
    std::cout << "Dielectric profile" << std::endl;
    for (unsigned int k=0; k<udisc.size()-1; k++)
      std::cout << "dielectric constant " << k+1 << " " << f[k] <<std::endl;
  }

  // Setting the value of D.

  Complex D = udisc.back() - udisc.front();
  if (extend)
    D *= 2.0;

  // Calculation of the fourier components.
  
  Complex K = 2.*pi/D; 
  cVector result(2*M+1,fortranArray);

  for (int m=0; m<=M; m++)
  {
    Complex E_m = 0.0;
      
    for (unsigned int k=0; k<int(disc.size()-1); k++)
    {
      const Complex d_k = udisc[k+1]-udisc[k];
      const Complex a   = (disc[k+1]-disc[k])/d_k;
   
      Complex factor;
      if (m == 0)
      {    
        factor = d_k*a;
        if (extend)
          factor *= 2.0;
      }
      else
      {
        const Complex t = -I*Real(m)*K;
        factor = (exp(t*udisc[k+1]) - exp(t*udisc[k])) / t;

        if (!extend)
        {
          Complex E = 2.*D/Real(m)/d_k;

          if (abs(1.0 - abs(E)) > 1e-12)
            factor *= a*(1.0 - eta/(1.0 - E*E));
          else
          {
            Complex A = (exp(2.* t*udisc[k+1])-exp(2.* t*udisc[k]))/(2.*t);
            factor = a*(factor- eta/2.*((udisc[k+1]-udisc[k])*exp(t*udisc[k])
               + A*exp(-t*udisc[k])));
          }
        }
        else
        {         
          factor += (exp(-t*udisc[k]) - exp(-t*udisc[k+1])) / t;

          Complex E, A;

          // Since the beginning and the end of the disc or udisc profile 
          // are not actually physical discontinuities (at the beginning the 
          // structure is mirrored, at the end the structure ends with an
          // electric or magnetic wall) the resolution should not be enhanced 
          // near these points. That is why the parametric function is 
          // slightly different for the first and the last segment.

          if ((k == 0) || (k == disc.size()-2))
            E = D/2.0/Real(m)/d_k;
          else
            E = D/Real(m)/d_k;

          if (abs(1.0 - abs(E)) > 1e-12)
            factor *= a*(1.0 - eta/(1.0 - E*E));
          else
          {
            if (k == 0)
            {
              A = (exp( 2.*t*udisc[k+1]) - exp( 2.*t*udisc[k]  )) / (2.*t) +
                  (exp(-2.*t*udisc[k])   - exp(-2.*t*udisc[k+1])) / (2.*t);

              factor = a*(factor-eta*((udisc[k+1]-udisc[k])*exp(-t*udisc[k+1])
                                       + A*exp(t*udisc[k+1])/2.));
            }
            else
            {
              A = (exp( 2.*t*udisc[k+1]) - exp( 2.*t*udisc[k]  )) / (2.*t) + 
                  (exp(-2.*t*udisc[k])   - exp(-2.*t*udisc[k+1])) / (2.*t);
              factor = a*(factor-eta*((udisc[k+1]-udisc[k])*exp(t*udisc[k])
                   + A*exp(-t*udisc[k])/2.));
            }
          }
        }
      }

      E_m += factor * f[k];
    }

    result(-m+M+1) = result(m+M+1) = E_m/D;

    // TODO_PDB: isn't this symmetry lost if extend == false?
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_ASR_pseudo
//
//  This is the 1D fourier integration in the new ASR basis, the difference 
//  with the previous definition is that this function takes u_disc as 
//  input, the length of the structure, and some booleans. The reason for 
//  this is that this function is only used for calcualting pseudo-fourier 
//  integrals and thus needs to know if we are at the beginning or the end 
//  of our stack.
//
//  TODO: join with previous function.
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier_ASR_pseudo(const vector<Complex>& f,
                           const vector<Complex>& disc, 
                           const vector<Complex>& udisc,
                           int M, const Complex* d, bool extend,
                           bool begin, bool end, const Real eta)
{
  if (0)
  {
    std::cout << "Fourier_ASR_pseudo" <<std::endl;
    std::cout << "Begin " << begin << ", end " << end <<std::endl;
    std::cout << "Input parameters, n_ASR " << eta << std::endl;
    for (int j=0; j<f.size(); j++)
      std::cout << "f[" << j << "] " << f[j] << std::endl;
    for (int j=0; j<disc.size(); j++)
      std::cout << "disc[" << j << "] " << disc[j] << std::endl;
    for (int j=0; j<udisc.size(); j++)
      std::cout << "udisc[" << j << "] " << udisc[j] << std::endl;
    std::cout << "d " << *d << std::endl;
    std::cout << "Boolean extend " << extend << std::endl;
  }

  // Setting the value of D.

  Complex D = *d;
 
  if (extend)
    D *= 2.0;

  // Calculation of the fourier components.
  
  Complex K = 2.*pi/D; 
  cVector result(2*M+1,fortranArray);

  for (int m=0; m<=M; m++)
  {
    Complex E_m = 0.0;
   
    const Complex d_k = udisc[1]-udisc[0];
    const Complex a   = (disc[1]-disc[0])/d_k;
   
    Complex factor;
    if (m == 0)
    {
      factor = d_k*a;
      if (extend)
        factor *= 2.0;
    }
    else
    {
      const Complex t = -I*Real(m)*K;
      factor = (exp(t*udisc[1]) - exp(t*udisc[0])) / t;

      if (! extend)
      {
        Complex E = 2.*D/Real(m)/d_k;

        if (abs(1.0 - abs(E)) > 1e-12)
          factor *= a*(1.0 - eta/(1.0 - E*E));
        else
        {
          Complex A = (exp(2.* t*udisc[1])-exp(2.* t*udisc[0]))/(2.*t);
          factor = a*(factor- eta/2.*((udisc[1]-udisc[0])*exp(t*udisc[0])
                                         + A*exp(-t*udisc[0])));
        }
      }
      else
      {
        factor += (exp(-t*udisc[0]) - exp(-t*udisc[1])) / t;
        
        Complex E, A;

        // Since the beginning and the end of the xdisc or udisc profile
        // are not actually physical discontinuities (at the beginning
        // the structure is mirrored, at the end the structure ends with
        // an electric or magnetic wall) the resolution should not
        // be enhanced near these points. That is why the parametric
        // function is slightly different for the first and the last
        // segment .

        if (begin || end)
          E = D/2.0/Real(m)/d_k;
        else
          E = D/Real(m)/d_k;

        if (abs(1.0 - abs(E)) > 1e-12)
          factor *= a*(1.0 - eta/(1.0 - E*E));
        else
        {
          if (begin)
          {
            A = (exp( 2.*t*udisc[1]) - exp( 2.*t*udisc[0])) / (2.*t) +
                (exp(-2.*t*udisc[0]) - exp(-2.*t*udisc[1])) / (2.*t);
            factor = a*(factor-eta*((udisc[1]-udisc[0])*exp(-t*udisc[1])
                                       + A*exp(t*udisc[1])/2.));
            }
          else
          {
            A = (exp( 2.*t*udisc[1]) - exp( 2.*t*udisc[0])) / (2.*t) +
                (exp(-2.*t*udisc[0]) - exp(-2.*t*udisc[1])) / (2.*t);
            factor = a*(factor-eta*((udisc[1]-udisc[0])*exp(t*udisc[0])
                                       + A*exp(-t*udisc[0])/2.));
          }
        }
      }
    }

    E_m += factor * f[0];

    result(-m+M+1) = result(m+M+1) = E_m/D;

    // TODO_PDB: isn't this symmetry lost if extend == false?
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D(const vector<Complex>& disc_x,
                   const vector<vector<Complex> >& disc_y,
                   const vector<vector<Complex> >& f,
                   int M, int N, bool extend)
{
  const Complex Lx = disc_x.back() - disc_x.front();

  cMatrix result(2*M+1,2*N+1,fortranArray);
  result = 0.0;
  
  for (int i=0; i<disc_x.size()-1; i++)
  {
    // Calculate 1D fourier transform of f(y) profile in this strip.

    cVector fourier_1D_y(2*N+1,fortranArray);   
    fourier_1D_y = fourier(f[i], disc_y[i], N, NULL, extend);

    // Calculate pseudo 1D fourier transform in x direction.

    vector<Complex> disc_i_x;
    disc_i_x.push_back(disc_x[i]);
    disc_i_x.push_back(disc_x[i+1]);

    vector<Complex> f_i_x; 
    f_i_x.push_back(1.0);

    cVector fourier_1D_x(2*M+1,fortranArray);
    fourier_1D_x = fourier(f_i_x, disc_i_x, M, &Lx, extend);

    // Calculate 2D fourier transform.

    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
        result(m+M+1,n+N+1) += fourier_1D_x(m+M+1) * fourier_1D_y(n+N+1);
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D_ASR
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_ASR(const vector<Complex>& disc_x,
                       const vector<vector<Complex> >& disc_y,
                       const vector<Complex>& disc_u,
                       const vector<vector<Complex> >& disc_v,
                       const vector<vector<Complex> >& f,
                       int M, int N, bool extend,
                       const Real eta_ASR)
{
  const Complex Lu = disc_u.back() - disc_u.front();

  cMatrix result(2*M+1,2*N+1,fortranArray);
  result = 0.0;
  
  for (int i=0; i<disc_u.size()-1; i++)
  {
    // Calculate 1D fourier transform of f(y) profile in this strip.

    cVector fourier_1D_v(2*N+1,fortranArray);
    fourier_1D_v = fourier_ASR(f[i], disc_y[i],disc_v[i], N, extend, eta_ASR);
    
    if (0)
    {
      std::cout << "Fourier_2D_ASR" << std::endl;
      std::cout << "Calculating the 1D fourier transform of f(y)" << std::endl;
      std::cout << "Input parameters, n_ASR " << eta_ASR << std::endl;
      for (int j=0; j<f[i].size(); j++)
        std::cout << "f[" << j << "] "<<f[i][j]<<std::endl;
      for (int j=0; j<disc_y[i].size(); j++)
        std::cout << "disc_y[" << j << "] " << disc_y[i][j] << std::endl;
      for (int j=0; j<disc_y[i].size(); j++)
        std::cout << "disc_v[" << j << "] " << disc_v[i][j] << std::endl;    
      for (int j=0; j<fourier_1D_v.size(); j++)
        std::cout << "fourier_1D_v[" << j << "] " << fourier_1D_v[j] 
                  << std::endl;
    }

    // Calculate pseudo 1D fourier transform in u direction.

    vector<Complex> disc_i_x;
    disc_i_x.push_back(disc_x[i]);
    disc_i_x.push_back(disc_x[i+1]);

    vector<Complex> disc_i_u;
    disc_i_u.push_back(disc_u[i]);
    disc_i_u.push_back(disc_u[i+1]);

    vector<Complex> f_i_x; 
    f_i_x.push_back(1.0);

    // Calculate booleans begin and end. These are necesary for the pseudo 
    // fourier ASR solver, because at the beginning and the end the parametric
    // formulation is slightly different, so the pseudo fourier solver should 
    // know when the first or the last piece will be calculated.
    
    bool begin = (i == 0);
    bool end   = (i == disc_x.size()-2);
  
    cVector fourier_1D_u(2*M+1,fortranArray);

    fourier_1D_u = fourier_ASR_pseudo(f_i_x, disc_i_x, disc_i_u, M, &Lu, 
                                      extend, begin, end, eta_ASR);
	
    if (0)
    {
      std::cout << "Calculating 1D transform in u direction" << std::endl;
      std::cout << "Input parameters, n_ASR " << eta_ASR << std::endl;
      std::cout << "f_i_x[" << 0 << "] " << f_i_x[0] << std::endl;
      std::cout << "disc_i_x[0] " << disc_i_x[0] << std::endl;
      std::cout << "disc_i_x[1] " << disc_i_x[1] << std::endl;
      std::cout << "disc_i_u[0] " << disc_i_u[0] << std::endl;
      std::cout << "disc_i_u[1] " << disc_i_u[1] << std::endl;
      
      std::cout << "Fourier 1D u" << std::endl;
      for (int j=0; j<fourier_1D_u.size(); j++)
        std::cout << "fourier_1D_u[" << j << "] " << fourier_1D_u[j]
                  << std::endl;
    }

    // Calculate 2D fourier transform.

    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
        result(m+M+1,n+N+1) += fourier_1D_u(m+M+1) * fourier_1D_v(n+N+1);
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D_split
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_split(const vector<Complex>& disc_x,
                         const vector<vector<Complex> >& disc_y,
                         const vector<vector<Complex> >& f,
                         int M, int N, bool extend)
{
  const Complex Lx = disc_x.back() - disc_x.front();

  const int MN = (2*M+1)*(2*N+1); 
  cMatrix result(MN,MN,fortranArray);
  result = 0.0;

  for (int i=0; i<disc_x.size()-1; i++)
  {
    // Calculate 1D fourier transform of f(y) profile in this strip.

    cVector fourier_1D_y(4*N+1,fortranArray);   
    fourier_1D_y = fourier(f[i], disc_y[i], 2*N, NULL, extend);

    // Invert Toeplitz matrix.

    cMatrix f_toep(2*N+1,2*N+1,fortranArray);
    for (int i1=-N; i1<=N; i1++)
      for (int i2=-N; i2<=N; i2++)
        f_toep(i1+N+1,i2+N+1) = fourier_1D_y(i1-i2 + 2*N+1);

    cMatrix inv_f_toep(2*N+1,2*N+1,fortranArray);
    inv_f_toep.reference(invert(f_toep)); // TODO: exploit Toeplitz.

    // Calculate pseudo 1D fourier transform in x direction.

    vector<Complex> disc_i_x;
    disc_i_x.push_back(disc_x[i]);
    disc_i_x.push_back(disc_x[i+1]);

    vector<Complex> f_i_x; 
    f_i_x.push_back(1.0);

    cVector fourier_1D_x(4*M+1,fortranArray);
    fourier_1D_x = fourier(f_i_x, disc_i_x, 2*M, &Lx, extend);

    // Fill result matrix.

    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
      {
        int i1 = (n+N+1) + (m+M)*(2*N+1);
      
        for (int j=-M; j<=M; j++)
          for (int l=-N; l<=N; l++)
          {
            int i2 = (l+N+1) + (j+M)*(2*N+1);

            result(i1,i2) += fourier_1D_x(m-j + 2*M+1) *
                               inv_f_toep(n   + N+1, l + N+1);
          }
      } 
  }

  return result;
}



/////////////////////////////////////////////////////////////////////////////
//
// fourier_2D_split_ASR
//
/////////////////////////////////////////////////////////////////////////////

cMatrix fourier_2D_split_ASR(const vector<Complex>& disc_x,
                         	const vector<vector<Complex> >& disc_y,
                         	const vector<Complex>& disc_u,
                         	const vector<vector<Complex> >& disc_v,
                         	const vector<vector<Complex> >& f,
                         	int M, int N, bool extend, 
                         	const Real eta)
{
  const Complex Lu = disc_u.back() - disc_u.front();
  
  const int MN = (2*M+1)*(2*N+1); 
  cMatrix result(MN,MN,fortranArray);
  result = 0.0;

  for (int i=0; i<disc_u.size()-1; i++)
  {
    // Calculate 1D fourier transform of f(v) profile in this strip.

    if (0)
    {      
      std::cout << "Calculating 1D fourier transform of f(y)" << std::endl;
      std::cout << "Input parameters, n_ASR " << eta << std::endl;
      for (int j=0; j<f[i].size(); j++)
        std::cout << "f[" << j << "] " << f[i][j] << std::endl;
      for (int j=0; j<disc_y[i].size(); j++)
        std::cout << "disc_y["<<j<<"] "<< disc_y[i][j] << std::endl;
      for(int j=0; j<disc_y[i].size(); j++)
        std::cout << "disc_v["<<j<<"] " << disc_v[i][j] << std::endl;
    }

    cVector fourier_1D_v(4*N+1,fortranArray);
	
    fourier_1D_v = fourier_ASR(f[i],disc_y[i],disc_v[i],2*N,extend,eta);
	
    // Invert Toeplitz matrix.

    cMatrix f_toep(2*N+1,2*N+1,fortranArray);
    for (int i1=-N; i1<=N; i1++)
      for (int i2=-N; i2<=N; i2++)
        f_toep(i1+N+1,i2+N+1) = fourier_1D_v(i1-i2 + 2*N+1);

    cMatrix inv_f_toep(2*N+1,2*N+1,fortranArray);
    inv_f_toep.reference(invert(f_toep)); // TODO: exploit Toeplitz.

    // Calculate pseudo 1D fourier transform in u direction.

    vector<Complex> disc_i_x;
    disc_i_x.push_back(disc_x[i]);
    disc_i_x.push_back(disc_x[i+1]);

    vector<Complex> disc_i_u;
    disc_i_u.push_back(disc_u[i]);
    disc_i_u.push_back(disc_u[i+1]);

    vector<Complex> f_i_x; 
    f_i_x.push_back(1.0);

    bool begin = (i == 0);
    bool end   = (i == disc_x.size()-2);
     
    cVector fourier_1D_u(4*M+1,fortranArray);
    
    if (0)
    {
      std::cout << "Calculating 1D transform in u direction" << std::endl;
      std::cout << "Input parameters, n_ASR " << eta << std::endl;
      std::cout << "f_i_x[" << 0 <<"] " << f_i_x[0] << std::endl;
      std::cout << "disc_i_x[0] "<< disc_i_x[0] << std::endl;
      std::cout << "disc_i_x[1] "<< disc_i_x[1] << std::endl;
      std::cout << "disc_i_u[0] "<< disc_i_u[0] << std::endl;
      std::cout << "disc_i_u[1] "<< disc_i_u[1] << std::endl;
    }

    fourier_1D_u = fourier_ASR_pseudo(f_i_x, disc_i_x, disc_i_u,
                                      2*M, &Lu, extend, begin, end, eta);
	
    // Fill result matrix.

    for (int m=-M; m<=M; m++)
      for (int n=-N; n<=N; n++)
      {
        int i1 = (n+N+1) + (m+M)*(2*N+1);
      
        for (int j=-M; j<=M; j++)
          for (int l=-N; l<=N; l++)
          {
            int i2 = (l+N+1) + (j+M)*(2*N+1);

            result(i1,i2) += fourier_1D_u(m-j + 2*M+1) *
                               inv_f_toep(n   + N+1, l + N+1);
          }
      } 
  }

  return result;
}
