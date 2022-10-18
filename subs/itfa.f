C*********************************************************************
C subroutine Inhomogeneous Temperature Function Add ******************
C            -             -           -        -   ******************
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
C For the coefficients CA, CB and CC (see ITFCF) ITFA will return    C
C derivatives 0 to IHD of f( r ) with f[nd]( r ) in DERV( nd + 1 )   C
C                                                                    C
C Currently IHD may be no larger than 4 although the routine could   C
C easily be modified for higher derivatives if the need arose.       C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IHD      : The number of the highest derivative requested.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD      : Current radius.                                     C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C     CA       : Coefficient of term in f(r). See above.             C
C     CB       : Coefficient of term in f(r). See above.             C
C     CC       : Coefficient of term in f(r). See above.             C
C     DERV     : Array of length atleast ( IHD + 1 ).                C
C                 DERV( nd + 1 ) returned with nd^{th} deriv. of f.  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IHD
      DOUBLE PRECISION RAD, RI, RO, CA, CB, CC, DERV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION PI, ROMRI, RMRI, LOW, FAC, TERMA, TERMB,
     1                 DOPRND
      PARAMETER ( PI=3.14159265358979312D0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF (      DABS( CA ).LT.LOW      .AND.
     1          DABS( CB ).LT.LOW      .AND.
     2          DABS( CC ).LT.LOW     )      RETURN
C     .
      RMRI  = RAD - RI
      ROMRI = RO - RI
      IF ( DABS( ROMRI ).LT.LOW ) THEN
        PRINT *,' Subroutine ITFA.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check on bounds for IHD
C     .
      IF ( IHD.LT.0 .OR. IHD.GT.4 ) THEN
        PRINT *,' Subroutine ITFA.'
        PRINT *,' IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC    = 0.5d0*PI/ROMRI
      DOPRND = FAC*RMRI
C
      TERMA  = CA*DSIN( DOPRND )
      TERMB  = CB*DCOS( DOPRND )*(-1.0d0)/FAC
      DERV( 1 ) = DERV( 1 ) + TERMA + TERMB + CC
      IF ( IHD.EQ.0 ) RETURN
C     .
      TERMA  = CA*FAC*DCOS( DOPRND )
      TERMB  = CB*DSIN( DOPRND )
      DERV( 2 ) = DERV( 2 ) + TERMA + TERMB
      IF ( IHD.EQ.1 ) RETURN
C     .
      TERMA  = (-1.0d0)*CA*FAC*FAC*DSIN( DOPRND )
      TERMB  = CB*FAC*DCOS( DOPRND )
      DERV( 3 ) = DERV( 3 ) + TERMA + TERMB
      IF ( IHD.EQ.2 ) RETURN
C     .
      TERMA  = (-1.0d0)*CA*FAC*FAC*FAC*DCOS( DOPRND )
      TERMB  = (-1.0d0)*CB*FAC*FAC*DSIN( DOPRND )
      DERV( 4 ) = DERV( 4 ) + TERMA + TERMB
      IF ( IHD.EQ.3 ) RETURN
C     .
      TERMA  = CA*FAC*FAC*FAC*FAC*DSIN( DOPRND )
      TERMB  = (-1.0d0)*CB*FAC*FAC*FAC*DCOS( DOPRND )
      DERV( 5 ) = DERV( 5 ) + TERMA + TERMB
      IF ( IHD.EQ.4 ) RETURN
C     .
      STOP
      END
C*********************************************************************

