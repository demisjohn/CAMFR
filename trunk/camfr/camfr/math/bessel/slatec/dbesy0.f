*DECK DBESY0
      DOUBLE PRECISION FUNCTION DBESY0 (X)
C***BEGIN PROLOGUE  DBESY0
C***PURPOSE  Compute the Bessel function of the second kind of order
C            zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (BESY0-S, DBESY0-D)
C***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESY0(X) calculates the double precision Bessel function of the
C second kind of order zero for double precision argument X.
C
C Series for BY0        on the interval  0.          to  1.60000E+01
C                                        with weighted error   8.14E-32
C                                         log weighted error  31.09
C                               significant figures required  30.31
C                                    decimal places required  31.73
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B0MP, DBESJ0, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESY0
      DOUBLE PRECISION X, BY0CS(19), AMPL, THETA, TWODPI, XSML,
     1  Y, D1MACH, DCSEVL, DBESJ0
      LOGICAL FIRST
      SAVE BY0CS, TWODPI, NTY0, XSML, FIRST
      DATA BY0CS(  1) / -.1127783939 2865573217 9398054602 8 D-1      /
      DATA BY0CS(  2) / -.1283452375 6042034604 8088453183 8 D+0      /
      DATA BY0CS(  3) / -.1043788479 9794249365 8176227661 8 D+0      /
      DATA BY0CS(  4) / +.2366274918 3969695409 2415926461 3 D-1      /
      DATA BY0CS(  5) / -.2090391647 7004862391 9622395034 2 D-2      /
      DATA BY0CS(  6) / +.1039754539 3905725209 9924657638 1 D-3      /
      DATA BY0CS(  7) / -.3369747162 4239720967 1877534503 7 D-5      /
      DATA BY0CS(  8) / +.7729384267 6706671585 2136721637 1 D-7      /
      DATA BY0CS(  9) / -.1324976772 6642595914 4347606896 4 D-8      /
      DATA BY0CS( 10) / +.1764823261 5404527921 0038936315 8 D-10     /
      DATA BY0CS( 11) / -.1881055071 5801962006 0282301206 9 D-12     /
      DATA BY0CS( 12) / +.1641865485 3661495027 9223718574 9 D-14     /
      DATA BY0CS( 13) / -.1195659438 6046060857 4599100672 0 D-16     /
      DATA BY0CS( 14) / +.7377296297 4401858424 9411242666 6 D-19     /
      DATA BY0CS( 15) / -.3906843476 7104373307 4090666666 6 D-21     /
      DATA BY0CS( 16) / +.1795503664 4361579498 2912000000 0 D-23     /
      DATA BY0CS( 17) / -.7229627125 4480104789 3333333333 3 D-26     /
      DATA BY0CS( 18) / +.2571727931 6351685973 3333333333 3 D-28     /
      DATA BY0CS( 19) / -.8141268814 1636949333 3333333333 3 D-31     /
      DATA TWODPI / 0.6366197723 6758134307 5535053490 057 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESY0
      IF (FIRST) THEN
         NTY0 = INITDS (BY0CS, 19, 0.1*REAL(D1MACH(3)))
         XSML = SQRT(4.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBESY0',
     +   'X IS ZERO OR NEGATIVE', 1, 2)
      IF (X.GT.4.0D0) GO TO 20
C
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBESY0 = TWODPI*LOG(0.5D0*X)*DBESJ0(X) + .375D0 + DCSEVL (
     1  .125D0*Y-1.D0, BY0CS, NTY0)
      RETURN
C
 20   CALL D9B0MP (X, AMPL, THETA)
      DBESY0 = AMPL * SIN(THETA)
      RETURN
C
      END
