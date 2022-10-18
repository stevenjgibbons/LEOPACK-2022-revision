C*********************************************************************
C subroutine Inhomogeneous Temperature Function Coefficient Find *****
C            -             -           -        -           -    *****
C Steve Gibbons Mon Jan 17 11:25:35 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C The function f( r ) for ri .le. r .le. ro is defined by            C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C Then clearly,    f( ri )  = CB 2 ( ri-ro )/pi + CC                 C
C                  f'( ri ) = CA 0.5 pi/(ro-ri)                      C
C                                                                    C
C                  f( ro )  = CA                + CC                 C
C                  f'( ro ) = CB                                     C
C                                                                    C
C ITFCF chooses the coefficients CA, CB and CC such that f will      C
C have the desired properties at ri and ro.                          C
C                                                                    C
C Set KIB (KOB) to 1 to fix the temperature at the inner (outer) bnd C
C Set KIB (KOB) to 2 to fix the heat flux at the inner (outer) bnd.  C
C                                                                    C
C VALIB (VALOB) is the actual value which f( r ) or f'( r ) is       C
C to attain at the inner (outer) boundary, depending upon the        C
C setting of KIB (KOB).                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KIB      : 1 if T( ri ) is to be set to the value VALIB.       C
C                2 if dT/dr( ri ) is to be set to the value VALIB.   C
C                                                                    C
C     KOB      : 1 if T( ro ) is to be set to the value VALOB.       C
C                2 if dT/dr( ro ) is to be set to the value VALOB.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C                                                                    C
C     VALIB    : Value for temp./temp. gradient at inner boundary.   C
C     VALOB    : Value for temp./temp. gradient at outer boundary.   C
C                                                                    C
C     CA       : Coefficient of term in f(r). See above.             C
C     CB       : Coefficient of term in f(r). See above.             C
C     CC       : Coefficient of term in f(r). See above.             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITFCF( KIB, KOB, RI, RO, VALIB, VALOB, CA, CB, CC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KIB, KOB
      DOUBLE PRECISION RI, RO, VALIB, VALOB, CA, CB, CC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION PI, ROMRI, LOW
      PARAMETER ( PI=3.14159265358979312D0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      ROMRI = RO - RI
      IF ( DABS( ROMRI ).LT.LOW ) THEN
        PRINT *,' Subroutine ITFCF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . First deal with the case of fixed temperature
C     . at both inner and outer boundaries
C     .
      IF ( KIB.EQ.1 .AND. KOB.EQ.1 ) THEN
        CA = VALOB
        CB = VALIB*PI*(-0.5d0)*ROMRI
        CC = 0.0d0
        RETURN
      ENDIF
C     .
C     . Fixed temperature at inner boundary
C     . Fixed flux at outer boundary
C     .
      IF ( KIB.EQ.1 .AND. KOB.EQ.2 ) THEN
        CA = 0.0d0
        CB = VALOB
        CC = VALIB - CB*2.0d0*(RI-RO)/PI
        RETURN
      ENDIF
C     .
C     . Fixed flux at inner boundary
C     . Fixed temperature at outer boundary
C     .
      IF ( KIB.EQ.2 .AND. KOB.EQ.1 ) THEN
        CA = VALIB*ROMRI*2.0d0/PI
        CB = 0.0d0
        CC = VALOB - CA
        RETURN
      ENDIF
C     .
C     . Fixed flux at inner and outer boundaries
C     .
      IF ( KIB.EQ.2 .AND. KOB.EQ.2 ) THEN
        CA = VALIB*ROMRI*2.0d0/PI
        CB = VALOB
        CC = 0.0d0
        RETURN
      ENDIF
C     .
      PRINT *,' Subroutine ITFCF.'
      PRINT *,' KIB = ', KIB,' KOB = ',KOB
      PRINT *,' Illegal values for KIB and KOB.'
      PRINT *,' Program aborted.'
C     .
      STOP
      END
C*********************************************************************

