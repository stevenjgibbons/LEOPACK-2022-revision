C*********************************************************************
C subroutine Poloidal Magnetic Field Alpha and Beta Find *************
C            -        -        -     -         -    -    *************
C Steve Gibbons Mon Oct 11 11:59:12 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If a poloidal radial function is expanded in the form             C
C                                                                    C
C  P( r ) = sum_i cos ( a_i r - b_i )                                C
C                                                                    C
C  when both inner and outer boundaries are insulating,              C
C  then for L, RI and RO (spherical harmonic degree, radius of       C
C  inner and outer boundary respectively) then a ( ALPHA ) and b     C
C  ( BETA ) must satisfy equations F1 = 0 and F2 = 0 as described    C
C  in Function PMFSEF.                                               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ALPHA     : Value a_i. Root to F1 and F2 in PMFSEF.            C
C                 On input is a guess for the i^{th} such root -     C
C                 This means that ALPHA is returned lying between    C
C                 i * pi/dist and ( i - 1 )*pi/dist.                 C
C                  ( dist = ro - ri )                                C
C                                                                    C
C     BETA      : Corresponding b_i value. Not constrained between   C
C                 any limits as is unique only up to periodicity.    C
C                                                                    C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     I         : Number of root. See above.                         C
C                                                                    C
C     L         : Spherical harmonic degree.                         C
C                                                                    C
C     ITMX      : Maximum iterations allowed to find a_i and b_i.    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PMFABF( ALPHA, BETA, I, RI, RO, L, ITMX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, L, ITMX
      DOUBLE PRECISION ALPHA, BETA, RI, RO
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NEQN
      PARAMETER ( NEQN = 2 )
      DOUBLE PRECISION PMFSEF, PI, ERR, XVEC( NEQN ), WORKV( NEQN ),
     1                 WORKA( NEQN, NEQN ), DPRARR( 2 ), PIN,
     2                 PINM1, DELPI, DIST
      INTEGER INFO, IWORK( NEQN ), INTARR( 1 ), NOIT
      EXTERNAL PMFSEF
      PARAMETER ( PI=3.14159265358979312D0, ERR = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( BETA.EQ.0.0d0 ) BETA = 1.5d0
C
      IF ( ITMX.LT.2 ) THEN
         PRINT *,' Subroutine PMFABF. ITMX = ', ITMX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Fill arrays for sending to MNEWTR
C
      DPRARR( 1 ) = RI
      DPRARR( 2 ) = RO
      DIST = RO - RI
C
      INTARR( 1 ) = L
C
      PIN    = DBLE( I )*PI/DIST
      PINM1  = (DBLE( I - 1 )*PI + ERR)/DIST
C
C ALPHA and BETA are our initial guesses -
C there is no constraint upon BETA, although ALPHA
C must lie between ( I - 1 )*PI and I*PI
C otherwise, we will ignore this guess and loop around
C lot's of random guesses - this will be at line 50.
C
      IF ( ALPHA.LE.PINM1 .OR. ALPHA.GT.PIN ) GOTO 50
C
C Ok so our guess for ALPHA seems o.k.
C Let's try and solve using this and our guess for BETA
C
      XVEC( 1 ) = ALPHA
      XVEC( 2 ) = BETA
C
      CALL MNEWTR( PMFSEF, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1             WORKA, IWORK, INTARR, DPRARR )
C
C Check to see if MNEWTR has found a (possibly a useless)
C solution ...
C
      IF ( INFO.EQ.-1 ) THEN
        PRINT *,' In MNEWTR, the LAPACK routine DGETRF'
        PRINT *,' returned a value of ', IWORK( 1 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INFO.EQ.-2 ) THEN
        PRINT *,' In MNEWTR, the LAPACK routine DGETRS'
        PRINT *,' returned a value of ', IWORK( 1 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INFO.EQ.-3 ) THEN
        PRINT *,' In MNEWTR, too many iterations were'
        PRINT *,' needed. ITMX = ', ITMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check our solution is in the right bounds
C If so, return ALPHA and BETA
C
      IF ( XVEC( 1 ).GT.PINM1 .AND. XVEC( 1 ).LE.PIN ) THEN
        ALPHA = XVEC( 1 )
        BETA  = XVEC( 2 )
        RETURN
      ENDIF
C
 50   CONTINUE
C
C O.k. our guesses weren't any good and so we will
C have to try and guess values ...
C
      DELPI = PI/(DBLE( ITMX + 1 )*DIST)
      DO NOIT = 1, ITMX
        XVEC( 1 ) = PIN - DBLE( NOIT )*DELPI
        XVEC( 2 ) = BETA
C
        CALL MNEWTR( PMFSEF, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1               WORKA, IWORK, INTARR, DPRARR )
C
C Check to see if MNEWTR has found a (possibly a useless)
C solution ...
C
        IF ( INFO.EQ.-1 ) THEN
          PRINT *,' In MNEWTR, the LAPACK routine DGETRF'
          PRINT *,' returned a value of ', IWORK( 1 )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        IF ( INFO.EQ.-2 ) THEN
          PRINT *,' In MNEWTR, the LAPACK routine DGETRS'
          PRINT *,' returned a value of ', IWORK( 1 )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        IF ( INFO.EQ.-3 ) THEN
          PRINT *,' In MNEWTR, too many iterations were'
          PRINT *,' needed. ITMX = ', ITMX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C o.k. - so A solution was found. Let's see
C if it is in range ...
C
        IF ( XVEC( 1 ).GT.PINM1 .AND. XVEC( 1 ).LE.PIN ) THEN
          ALPHA = XVEC( 1 )
          BETA  = XVEC( 2 )
          RETURN
        ENDIF
C
 60   CONTINUE
      ENDDO
C
      PRINT *,' All you guesses for ALPHA are '
      PRINT *,' now exhausted and your solution '
      PRINT *,' still cannot be found.'
      STOP
      END
C*********************************************************************
