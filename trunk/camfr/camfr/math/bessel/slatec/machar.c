/*

This file combines the single and double precision versions of machar,
selected by cc -DSP or cc -DDP.  This feature provided by D. G. Hough,
August 3, 1988.

Minor modifications by Peter Bienstman, to allow creating both single
and double precision versions at the same time.

*/

#undef REAL
#undef ZERO
#undef ONE
#undef PREC
#undef REALSIZE
#undef RMACHAR
#undef MACHAR

#ifdef SP

#define REAL float
#define ZERO 0.0
#define ONE 1.0
#define PREC "Single "
#define REALSIZE 1
#define RMACHAR rmachars

#ifdef FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#define MACHAR machars
#endif
#ifdef FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#define MACHAR machars_
#endif
#ifdef FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES
#define MACHAR machars__
#endif

#ifndef MACHAR
#define MACHAR machars_
#endif

#endif

#ifdef DP

#define REAL double
#define ZERO 0.0e0
#define ONE 1.0e0
#define PREC "Double "
#define REALSIZE 2
#define RMACHAR rmachard

#ifdef FORTRAN_SYMBOLS_WITHOUT_TRAILING_UNDERSCORES
#define MACHAR machard
#endif
#ifdef FORTRAN_SYMBOLS_WITH_SINGLE_TRAILING_UNDERSCORE
#define MACHAR machard_
#endif
#ifdef FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES
#define MACHAR machard__
#endif

#ifndef MACHAR
#define MACHAR machard_
#endif

#endif

#include <math.h>

#define ABS(xxx) ((xxx>ZERO)?(xxx):(-xxx))

void
RMACHAR (ibeta,it,irnd,ngrd,machep,negep,iexp,
         minexp,maxexp,eps,epsneg,xmin,xmax)

      int *ibeta,*iexp,*irnd,*it,*machep,*maxexp,*minexp,*negep,*ngrd;
      REAL *eps,*epsneg,*xmax,*xmin;

/*

    This subroutine is intended to determine the parameters of the
    floating-point arithmetic system specified below.  The
    determination of the first three uses an extension of an algorithm
    due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
    but not all, of the improvements suggested by M. Gentleman and S.
    Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
    program was published in the book Software Manual for the
    Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
    Englewood Cliffs, NJ, 1980.  The present program is a
    translation of the Fortran 77 program in W. J. Cody, "MACHAR:
    A subroutine to dynamically determine machine parameters".
    TOMS (14), 1988.

   Parameter values reported are as follows:

        ibeta   - the radix for the floating-point representation
        it      - the number of base ibeta digits in the floating-point
                  significand
        irnd    - 0 if floating-point addition chops
                  1 if floating-point addition rounds, but not in the
                    IEEE style
                  2 if floating-point addition rounds in the IEEE style
                  3 if floating-point addition chops, and there is
                    partial underflow
                  4 if floating-point addition rounds, but not in the
                    IEEE style, and there is partial underflow
                  5 if floating-point addition rounds in the IEEE style,
                    and there is partial underflow
        ngrd    - the number of guard digits for multiplication with
                  truncating arithmetic.  It is
                  0 if floating-point arithmetic rounds, or if it
                    truncates and only  it  base  ibeta digits
                    participate in the post-normalization shift of the
                    floating-point significand in multiplication;
                  1 if floating-point arithmetic truncates and more
                    than  it  base  ibeta  digits participate in the
                    post-normalization shift of the floating-point
                    significand in multiplication.
        machep  - the largest negative integer such that
                  1.0+FLOAT(ibeta)**machep .NE. 1.0, except that
                  machep is bounded below by  -(it+3)
        negeps  - the largest negative integer such that
                  1.0-FLOAT(ibeta)**negeps .NE. 1.0, except that
                  negeps is bounded below by  -(it+3)
        iexp    - the number of bits (decimal places if ibeta = 10)
                  reserved for the representation of the exponent
                  (including the bias or sign) of a floating-point
                  number
        minexp  - the largest in magnitude negative integer such that
                  FLOAT(ibeta)**minexp is positive and normalized
        maxexp  - the smallest positive power of  BETA  that overflows
        eps     - the smallest positive floating-point number such
                  that  1.0+eps .NE. 1.0. In particular, if either
                  ibeta = 2  or  IRND = 0, eps = FLOAT(ibeta)**machep.
                  Otherwise,  eps = (FLOAT(ibeta)**machep)/2
        epsneg  - A small positive floating-point number such that
                  1.0-epsneg .NE. 1.0. In particular, if ibeta = 2
                  or  IRND = 0, epsneg = FLOAT(ibeta)**negeps.
                  Otherwise,  epsneg = (ibeta**negeps)/2.  Because
                  negeps is bounded below by -(it+3), epsneg may not
                  be the smallest number that can alter 1.0 by
                  subtraction.
        xmin    - the smallest non-vanishing normalized floating-point
                  power of the radix, i.e.,  xmin = FLOAT(ibeta)**minexp
        xmax    - the largest finite floating-point number.  In
                  particular  xmax = (1.0-epsneg)*FLOAT(ibeta)**maxexp
                  Note - on some machines  xmax  will be only the
                  second, or perhaps third, largest number, being
                  too small by 1 or 2 units in the last digit of
                  the significand.

      Latest revision - August 4, 1988

      Author - W. J. Cody
               Argonne National Laboratory

*/

{
      int i,iz,j,k;
      int mx,itmp,nxres;
      REAL a,b,beta,betain,one,y,z,zero;
      REAL betah,t,tmp,tmpa,tmp1,two;

      (*irnd) = 1;
      one = (REAL)(*irnd);
      two = one + one;
      a = two;
      b = a;
      zero = 0.0e0;

/*
  Determine ibeta,beta ala malcolm.
*/

      tmp = ((a+one)-a)-one;

      while (tmp == zero) {
         a = a+a;
         tmp = a+one;
         tmp1 = tmp-a;
         tmp = tmp1-one;
      }

      tmp = a+b;
      itmp = (int)(tmp-a);
      while (itmp == 0) {
         b = b+b;
         tmp = a+b;
         itmp = (int)(tmp-a);
      }

      *ibeta = itmp;
      beta = (REAL)(*ibeta);

/*
  Determine irnd, it.
*/

      (*it) = 0;
      b = one;
      tmp = ((b+one)-b)-one;

      while (tmp == zero) {
         *it = *it+1;
         b = b*beta;
         tmp = b+one;
         tmp1 = tmp-b;
         tmp = tmp1-one;
      }

      *irnd = 0;
      betah = beta/two;
      tmp = a+betah;
      tmp1 = tmp-a;
      if (tmp1 != zero) *irnd = 1;
      tmpa = a+beta;
      tmp = tmpa+betah;
      if ((*irnd == 0) && (tmp-tmpa != zero)) *irnd = 2;

/*
  Determine negep, epsneg.
*/

      (*negep) = (*it) + 3;
      betain = one / beta;
      a = one;

      for (i = 1; i<=(*negep); i++) {
         a = a * betain;
      }

      b = a;
      tmp = (one-a);
      tmp = tmp-one;

      while (tmp == zero) {
         a = a*beta;
         *negep = *negep-1;
         tmp1 = one-a;
         tmp = tmp1-one;
      }

      (*negep) = -(*negep);
      (*epsneg) = a;

/*
  Determine machep, eps.
*/

      (*machep) = -(*it) - 3;
      a = b;
      tmp = one+a;

      while (tmp-one == zero) {
         a = a*beta;
         *machep = *machep+1;
         tmp = one+a;
      }

      *eps = a;

/*
  Determine ngrd.
*/

      (*ngrd) = 0;
      tmp = one+*eps;
      tmp = tmp*one;
      if (((*irnd) == 0) && (tmp-one) != zero) (*ngrd) = 1;

/*
  Determine iexp, minexp, xmin.

  Loop to determine largest i such that
         (1/beta) ** (2**(i))
  does not underflow.
  Exit from loop is signaled by an underflow.
*/

      i = 0;
      k = 1;
      z = betain;
      t = one+*eps;
      nxres = 0;

      for (;;) {
         y = z;
         z = y * y;

/*
  Check for underflow.
*/

         a = z * one;
         tmp = z*t;
         if ((a+a == zero) || (ABS(z) > y)) break;
         tmp1 = tmp*betain;
         if (tmp1*beta == z) break;
         i = i + 1;
         k = k+k;
      }

/*
  Determine k such that (1/beta)**k does not underflow.
  first set  k = 2 ** i.
*/

      (*iexp) = i + 1;
      mx = k + k;
      if (*ibeta == 10) {

/*
  For decimal machines only
*/

         (*iexp) = 2;
         iz = *ibeta;
         while (k >= iz) {
            iz = iz * (*ibeta);
            (*iexp) = (*iexp) + 1;
         }
         mx = iz + iz - 1;
      }

/*
  Loop to determine minexp, xmin.
  Exit from loop is signaled by an underflow.
*/

      for (;;) {
         (*xmin) = y;
         y = y * betain;
         a = y * one;
         tmp = y*t;
         tmp1 = a+a;
         if ((tmp1 == zero) || (ABS(y) >= (*xmin))) break;
         k = k + 1;
         tmp1 = tmp*betain;
         tmp1 = tmp1*beta;

         if ((tmp1 == y) && (tmp != y)) {
            nxres = 3;
            *xmin = y;
            break;
         }

      }

      (*minexp) = -k;

/*
  Determine maxexp, xmax.
*/

      if ((mx <= k+k-3) && ((*ibeta) != 10)) {
         mx = mx + mx;
         (*iexp) = (*iexp) + 1;
      }

      (*maxexp) = mx + (*minexp);

/*
  Adjust *irnd to reflect partial underflow.
*/

      (*irnd) = (*irnd)+nxres;

/*
  Adjust for IEEE style machines.
*/

      if ((*irnd) >= 2) (*maxexp) = (*maxexp)-2;

/*
  Adjust for machines with implicit leading bit in binary
  significand and machines with radix point at extreme
  right of significand.
*/

      i = (*maxexp) + (*minexp);
      if (((*ibeta) == 2) && (i == 0)) (*maxexp) = (*maxexp) - 1;
      if (i > 20) (*maxexp) = (*maxexp) - 1;
      if (a != y) (*maxexp) = (*maxexp) - 2;
      (*xmax) = one - (*epsneg);
      tmp = (*xmax)*one;
      if (tmp != (*xmax)) (*xmax) = one - beta * (*epsneg);
      (*xmax) = (*xmax) / (beta * beta * beta * (*xmin));
      i = (*maxexp) + (*minexp) + 3;
      if (i > 0) {

         for (j = 1; j<=i; j++ ) {
             if ((*ibeta) == 2) (*xmax) = (*xmax) + (*xmax);
             if ((*ibeta) != 2) (*xmax) = (*xmax) * beta;
         }

      }

    return;

}

void
MACHAR (REAL *xmin, REAL *xmax, REAL *epsneg, REAL *eps, REAL *log10_ibeta)
{ 
  int ibeta, iexp, irnd, it, machep, maxexp, minexp, negep, ngrd;

  RMACHAR (&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp,
	   &maxexp, eps, epsneg, xmin, xmax);

  *log10_ibeta = log10 ((REAL) ibeta);
}
