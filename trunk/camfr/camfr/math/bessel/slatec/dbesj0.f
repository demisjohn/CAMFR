*DECK DBESJ0
      DOUBLE PRECISION FUNCTION DBESJ0 (X)
C***BEGIN PROLOGUE  DBESJ0
C***PURPOSE  Compute the Bessel function of the first kind of order
C            zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (BESJ0-S, DBESJ0-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESJ0(X) calculates the double precision Bessel function of
C the first kind of order zero for double precision argument X.
C
C Series for BJ0        on the interval  0.          to  1.60000E+01
C                                        with weighted error   4.39E-32
C                                         log weighted error  31.36
C                               significant figures required  31.21
C                                    decimal places required  32.00
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B0MP, DCSEVL, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DBESJ0
      DOUBLE PRECISION X, BJ0CS(19), AMPL, THETA, XSML, Y, D1MACH,
     1  DCSEVL
      LOGICAL FIRST
      SAVE BJ0CS, NTJ0, XSML, FIRST
      DATA BJ0CS(  1) / +.1002541619 6893913701 0731272640 74 D+0     /
      DATA BJ0CS(  2) / -.6652230077 6440513177 6787578311 24 D+0     /
      DATA BJ0CS(  3) / +.2489837034 9828131370 4604687266 80 D+0     /
      DATA BJ0CS(  4) / -.3325272317 0035769653 8843415038 54 D-1     /
      DATA BJ0CS(  5) / +.2311417930 4694015462 9049241177 29 D-2     /
      DATA BJ0CS(  6) / -.9911277419 9508092339 0485193365 49 D-4     /
      DATA BJ0CS(  7) / +.2891670864 3998808884 7339037470 78 D-5     /
      DATA BJ0CS(  8) / -.6121085866 3032635057 8184074815 16 D-7     /
      DATA BJ0CS(  9) / +.9838650793 8567841324 7687486364 15 D-9     /
      DATA BJ0CS( 10) / -.1242355159 7301765145 5158970068 36 D-10    /
      DATA BJ0CS( 11) / +.1265433630 2559045797 9158272103 63 D-12    /
      DATA BJ0CS( 12) / -.1061945649 5287244546 9148175129 59 D-14    /
      DATA BJ0CS( 13) / +.7470621075 8024567437 0989155840 00 D-17    /
      DATA BJ0CS( 14) / -.4469703227 4412780547 6270079999 99 D-19    /
      DATA BJ0CS( 15) / +.2302428158 4337436200 5230933333 33 D-21    /
      DATA BJ0CS( 16) / -.1031914479 4166698148 5226666666 66 D-23    /
      DATA BJ0CS( 17) / +.4060817827 4873322700 8000000000 00 D-26    /
      DATA BJ0CS( 18) / -.1414383600 5240913919 9999999999 99 D-28    /
      DATA BJ0CS( 19) / +.4391090549 6698880000 0000000000 00 D-31    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESJ0
      IF (FIRST) THEN
         NTJ0 = INITDS (BJ0CS, 19, 0.1*REAL(D1MACH(3)))
         XSML = SQRT(8.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0D0) GO TO 20
C
      DBESJ0 = 1.0D0
      IF (Y.GT.XSML) DBESJ0 = DCSEVL (.125D0*Y*Y-1.D0, BJ0CS, NTJ0)
      RETURN
C
 20   CALL D9B0MP (Y, AMPL, THETA)
      DBESJ0 = AMPL * COS(THETA)
C
      RETURN
      END
