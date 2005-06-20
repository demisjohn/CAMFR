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
//  TODO_PDB: see if it's better to create udisc outside of this function.
//
/////////////////////////////////////////////////////////////////////////////

cVector fourier_ASR(const vector<Complex>& f,
                    const vector<Complex>& disc, int M, const Complex* d,
                    bool extend, bool stretching, const Real eta)
{  
  // Creation of the stretched refractive index profile vector udisc.

  

  vector<Complex> udisc;
  if (stretching == true)
  {  
    Complex u_begin = real(disc[1]-disc[0]);
    Complex u_end   = real(disc[disc.size()-1]-disc[disc.size()-2]);
    udisc.push_back(0.0);

    Complex u_step = (real(u_begin) >= real(u_end)) ? u_begin : u_end;
    for (unsigned int i=1; i<disc.size(); i++)
      udisc.push_back(Real(i)*u_step);

    // Adding the imaginary (PML) part of the disc vector to udisc.

    for (unsigned int i=1; i<disc.size(); i++)
      udisc[i] += imag(disc[i])*I;
  }
  else
  {
    udisc = disc;
  }

  // Setting the value of D.

  Complex D = udisc.back() - udisc.front();
  if (extend)
    D *= 2.0;

  //std::cout << "Original discontinuity profile" << std::endl;
  //for (unsigned int k=0; k<disc.size(); k++)
  //  std::cout << "discontinuity "  << k+1 << " " << disc[k] << std::endl;

  //std::cout << "New discontinuity profile" << std::endl;
  //for (unsigned int k=0; k<udisc.size(); k++)
  //  std::cout << "discontinuity " <<< k+1 << " " << udisc[k] <<std::endl;
  
  //std::cout << "Dielectric profile" << std::endl;
  //for (unsigned int k=0; k<udisc.size()-1; k++)
  //  std::cout << "dielectric constant " << k+1 << " " << f[k] <<std::endl;

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
      if (m==0)
      {    
        factor = d_k*a;
        if(extend)
          factor *= 2.0;
      }
      else
      {
        const Complex t = -I*Real(m)*K;
        factor = (exp(t*udisc[k+1]) - exp(t*udisc[k])) / t;

        if (! extend)
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

          // TODO_PDB: add comment.

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

    result(-m+M+1)=result(m+M+1) = E_m/D;


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

#if 0
  // Check result

  Complex Ly = slabs[0]->get_width();

  std::ofstream file("fourier.txt");
  for (Real x=-real(disc_x.back()); x<=real(disc_x.back()); x+=0.01)
    for (Real y=-real(Ly); y<=real(Ly); y+=0.01)
    {
      Complex t=0.0;
      for (int m=-M; m<=M; m+=1)
        for (int n=-N; n<=N; n+=1)
          t += result(m+M+1,n+N+1)
            *exp(I*Real(m)*2.*pi*x/Lx/2.)*exp(I*Real(n)*2.*pi*y/Ly/2.);
      file << x << " " << y << " " << real(t) << std::endl;
    }
  file.close();
#endif
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

    //std::cout << "new " << i << " " << fourier_1D_y << std::endl;

    // Invert Toeplitz matrix.

    cMatrix f_toep(2*N+1,2*N+1,fortranArray);
    for (int i1=-N; i1<=N; i1++)
      for (int i2=-N; i2<=N; i2++)
        f_toep(i1+N+1,i2+N+1) = fourier_1D_y(i1-i2 + 2*N+1);

    cMatrix inv_f_toep(2*N+1,2*N+1,fortranArray);
    inv_f_toep.reference(invert(f_toep)); // TODO: exploit Toeplitz.

    //std::cout << "new " << i << " " << inv_f_toep << std::endl;

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


