*DECK DBESJ1
      DOUBLE PRECISION FUNCTION DBESJ1 (X)
C***BEGIN PROLOGUE  DBESJ1
C***PURPOSE  Compute the Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (BESJ1-S, DBESJ1-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESJ1(X) calculates the double precision Bessel function of the
C first kind of order one for double precision argument X.
C
C Series for BJ1        on the interval  0.          to  1.60000E+01
C                                        with weighted error   1.16E-33
C                                         log weighted error  32.93
C                               significant figures required  32.36
C                                    decimal places required  33.57
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B1MP, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   910401  Corrected error in code which caused values to have the
C           wrong sign for arguments less than 4.0.  (WRB)
C***END PROLOGUE  DBESJ1
      DOUBLE PRECISION X, BJ1CS(19), AMPL, THETA, XSML, XMIN, Y,
     1  D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BJ1CS, NTJ1, XSML, XMIN, FIRST
      DATA BJ1CS(  1) / -.1172614151 3332786560 6240574524 003 D+0    /
      DATA BJ1CS(  2) / -.2536152183 0790639562 3030884554 698 D+0    /
      DATA BJ1CS(  3) / +.5012708098 4469568505 3656363203 743 D-1    /
      DATA BJ1CS(  4) / -.4631514809 6250819184 2619728789 772 D-2    /
      DATA BJ1CS(  5) / +.2479962294 1591402453 9124064592 364 D-3    /
      DATA BJ1CS(  6) / -.8678948686 2788258452 1246435176 416 D-5    /
      DATA BJ1CS(  7) / +.2142939171 4379369150 2766250991 292 D-6    /
      DATA BJ1CS(  8) / -.3936093079 1831797922 9322764073 061 D-8    /
      DATA BJ1CS(  9) / +.5591182317 9468800401 8248059864 032 D-10   /
      DATA BJ1CS( 10) / -.6327616404 6613930247 7695274014 880 D-12   /
      DATA BJ1CS( 11) / +.5840991610 8572470032 6945563268 266 D-14   /
      DATA BJ1CS( 12) / -.4482533818 7012581903 9135059199 999 D-16   /
      DATA BJ1CS( 13) / +.2905384492 6250246630 6018688000 000 D-18   /
      DATA BJ1CS( 14) / -.1611732197 8414416541 2118186666 666 D-20   /
      DATA BJ1CS( 15) / +.7739478819 3927463729 8346666666 666 D-23   /
      DATA BJ1CS( 16) / -.3248693782 1119984114 3466666666 666 D-25   /
      DATA BJ1CS( 17) / +.1202237677 2274102272 0000000000 000 D-27   /
      DATA BJ1CS( 18) / -.3952012212 6513493333 3333333333 333 D-30   /
      DATA BJ1CS( 19) / +.1161678082 2664533333 3333333333 333 D-32   /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESJ1
      IF (FIRST) THEN
         NTJ1 = INITDS (BJ1CS, 19, 0.1*REAL(D1MACH(3)))
C
         XSML = SQRT(8.0D0*D1MACH(3))
         XMIN = 2.0D0*D1MACH(1)
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0D0) GO TO 20
C
      DBESJ1 = 0.0D0
      IF (Y.EQ.0.0D0) RETURN
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBESJ1',
     +   'ABS(X) SO SMALL J1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBESJ1 = 0.5D0*X
      IF (Y.GT.XSML) DBESJ1 = X*(.25D0 + DCSEVL (.125D0*Y*Y-1.D0,
     1  BJ1CS, NTJ1) )
      RETURN
C
 20   CALL D9B1MP (Y, AMPL, THETA)
      DBESJ1 = SIGN (AMPL, X) * COS(THETA)
C
      RETURN
      END
