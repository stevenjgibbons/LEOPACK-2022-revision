C*********************************************************************
C subroutine Inhomogeneous Temperature Harmonic Coefficient Allocate *
C            -             -           -        -           -        *
C Steve Gibbons Wed Jan 19 08:46:34 GMT 2000              Routine    C
C____________________________________________________________________C
C                                                                    C
C ITHCAR receives arrays defining sets of harmonics; these are       C
C                                                                    C
C  MHT( ih ) = 1 for a poloidal velocity harmonic                    C
C  MHT( ih ) = 2 for a toroidal velocity harmonic                    C
C  MHT( ih ) = 3 for a temperature harmonic                          C
C  MHT( ih ) = 4 for a poloidal magnetic field harmonic              C
C  MHT( ih ) = 5 for a toroidal magnetic field harmonic              C
C                                                                    C
C  MHL( ih ) = l, degree of spherical harmonic                       C
C  MHM( ih ) = m (order) for cos ( m phi ) dependence or             C
C              -m for sin ( m phi ) dependence.                      C
C                                                                    C
C  A fourth array MHI is assigned values by ITHCAR: If MHT( ih ) = 3 C
C then ITHCAR will look at the arrays HMIB and HMOB (harmonic map    C
C for inner/outer boundary) and the flags KIB and KOB and decide     C
C what value of temperature/temp. gradient that harmonic radial      C
C function should achieve at that boundary.                          C
C                                                                    C
C MHIB and MHOB are indexed by INDSHC.                               C
C                                                                    C
C These values are supplied to the routine ITFCF as VALIB and VALOB  C
C and along with RI and RO will be used to provide the coefficients  C
C CA, CB and CC which define the inhomogeneous temperature function  C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C The coefficients CA, CB and CC are stored in the array             C
C CAFIT (coefficient array for inhomogeneous temperature) which has  C
C dimensions ( 3, NITHMX ) where NITHMX limits the possible number   C
C of temperature harmonics.                                          C
C                                                                    C
C NITH is the number of inhomogeneous temperature harmonics which    C
C have already been assigned coefficients by ITHCAR.                 C
C                                                                    C
C For example if IH is a harmonic radial function with MHT( IH ) = 3 C
C and ITFCF calculates that f(r) should take the coefficients        C
C CA, CB and CC; then                                                C
C                                                                    C
C  CA is stored in CAFIT( 1, IITH )                                  C
C  CB is stored in CAFIT( 2, IITH )                                  C
C  CC is stored in CAFIT( 3, IITH )                                  C
C                                                                    C
C and the index IITH is stored in MHI( IH )                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KIB      : 1 if T( r ) is to be fixed at the inner boundary.   C
C                2 if dT/dr(r) is to be fixed at the inner boundary. C
C                                                                    C
C     KOB      : 1 if T( r ) is to be fixed at the outer boundary.   C
C                2 if dT/dr(r) is to be fixed at the outer boundary. C
C                                                                    C
C     NITH     : Number of inhomogeneous temperature harmonics       C
C                 with coeff.s stored in CAFIT.                      C
C                                                                    C
C     NITHMX   : Limit on NITH.                                      C
C                                                                    C
C     NH       : Number of harmonics                                 C
C     MHT      : Dim( * ). Harmonic type - see above.                C
C     MHL      : Dim( * ). Harmonic degree - see above.              C
C     MHM      : Dim( * ). Harmonic order  - see above.              C
C     MHI      : Dim( * ). 2nd index of CAFIT where coeffs are storedC
C                                                                    C
C     LH       : Maximum spherical harmonic degree.                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C                                                                    C
C     HMIB     : Harmonic map for inner boundary. Dim( LH*(LH+2) )   C
C     HMIBMT   : Harmonic map for inner boundary monopole term       C
C     HMOB     : Harmonic map for outer boundary. Dim( LH*(LH+2) )   C
C     HMOBMT   : Harmonic map for outer boundary monopole term       C
C                                                                    C
C                 MHIB and MHOB can refer to either temperature      C
C                 or temperature gradient depending upon the         C
C                 values of KIB and KOB.                             C
C                                                                    C
C                 The arrays are indexed by INDSHC.                  C
C                                                                    C
C     CAFIT    : Coeff arr. for inhomog. temp. Dim( 3, NITHMX ).     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITHCAR( KIB, KOB, NITH, NITHMX, NH, MHT, MHL, MHM,
     1                   MHI, LH, RI, RO, HMIB, HMIBMT, HMOB, HMOBMT,
     2                   CAFIT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KIB, KOB, NITH, NITHMX, NH, MHT( * ), MHL( * ),
     1        MHM( * ), MHI( * ), LH
      DOUBLE PRECISION RI, RO, HMIB( LH*(LH+2) ), HMIBMT,
     1                  HMOB( LH*(LH+2) ), HMOBMT, CAFIT( 3, NITHMX )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHARM, L, M, ICS, IH, INDSHC, NIT
      DOUBLE PRECISION VALIB, VALOB, CA, CB, CC, DLOW
      LOGICAL OK
      PARAMETER ( DLOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF ( NITH.LT.0 ) THEN
        PRINT *,' Subroutine ITHCAR.'
        PRINT *,' NITH = ', NITH,' on entry.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      OK = .TRUE.
C     .
      DO IH = 1, NH
        IF ( MHT( IH ).NE.3 ) THEN
          MHI( IH ) = -1
          GOTO 50
        ENDIF
C       .
        L = MHL( IH )
        IF ( MHM( IH ).LT.0 ) THEN
          M   = -MHM( IH )
          ICS = 2
        ELSE
          M   = MHM( IH )
          ICS = 1
        ENDIF
        IHARM = INDSHC( L, M, ICS )
C       .
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine ITHCAR.'
          PRINT *,' Harmonic ',IH,' has L = ', L
          PRINT *,' Maximum degree = ', LH
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( IHARM.EQ.0 ) THEN
          VALIB = HMIBMT
          VALOB = HMOBMT
        ELSE
          VALIB = HMIB( IHARM )
          VALOB = HMOB( IHARM )
        ENDIF
C       .
        CALL ITFCF( KIB, KOB, RI, RO, VALIB, VALOB, CA, CB, CC )
C       .
C       . OK - we have calculated values
C       . for CA, CB and CC, so let's loop around
C       . CAFIT to see if we have already stored
C       . such values.
C       .
        DO NIT = 1, NITH
          IF ( DABS( CA - CAFIT( 1, NIT ) ).LT.DLOW .AND.
     1         DABS( CB - CAFIT( 2, NIT ) ).LT.DLOW .AND.
     2         DABS( CC - CAFIT( 3, NIT ) ).LT.DLOW   ) THEN
            MHI( IH ) = NIT
            GOTO 50
          ENDIF
        ENDDO
C       .
        CALL CNTRIC( NITH, NITHMX, OK )
C       .
        IF ( OK ) THEN
          CAFIT( 1, NITH ) = CA
          CAFIT( 2, NITH ) = CB
          CAFIT( 3, NITH ) = CC
          MHI( IH ) = NITH
        ENDIF
C       .
 50     CONTINUE
      ENDDO
C     .
      IF ( OK ) RETURN
      PRINT *,' Subroutine ITHCAR.'
      PRINT *,' Number of requested coefficients = ', NITH
      PRINT *,' NITHMX = ', NITHMX
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
