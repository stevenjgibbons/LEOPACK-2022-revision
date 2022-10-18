C*********************************************************************
C subroutine VECtor * 2 DiSplay **************************************
C            ---      - - -     **************************************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C Outputs an eye-readable breakdown of two solution vectors with     C
C the radial function of each component being headed by its          C
C description.                                                       C
C Note that they must have identical dimensions and harmonic         C
C descriptions!!                                                     C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDIM      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     LU        : File number for output                             C
C     IMODE     : Mode of choice for output                          C
C   IMODE = 1 lists every harmonic with full title and radial value  C
C                                                                    C
C   IZF         : = 1 to display ALL vectors                         C
C                 = 2 to display only non-zero harmonics             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC1     : DP vector of dimension ( NDIM )                    C
C     RVEC2     : DP vector of dimension ( NDIM )                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VEC2DS ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1                    RVEC1, RVEC2, LU, IMODE, IZF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH, LU, IMODE, IZF,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, RVEC1( NDIM ), RVEC2( NDIM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, RAD, HTOT, TOL
      INTEGER IRN, IH, IND
      CHARACTER *(11) CHPOL, CHTOR, CHTHE, CH11
      CHARACTER *(3) CHCOS, CHSIN, CH3
      PARAMETER ( TOL = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IZF.EQ.2 ) THEN
        WRITE ( LU, * ) 'Non-zero radial functions only,'
      ENDIF
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine VEC2DS. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
      CHPOL = 'Poloidal   '
      CHTOR = 'Toroidal   '
      CHTHE = 'Temperature'
      CHCOS = 'Cos'
      CHSIN = 'Sin'
      H = ( RO - RI ) / DBLE( NR - 1 )
C
      DO IH = 1, NH
        IF ( IZF.EQ.2 ) THEN
           HTOT = 0.0d0
           DO IRN = 1, NR
             IND = ( IRN - 1 )*NH + IH
             HTOT = HTOT + ABS( RVEC1( IND ) ) + ABS( RVEC2( IND ) )
           ENDDO
           IF ( HTOT.LT.TOL ) GOTO 100
        ENDIF
        IF ( MHT( IH ).EQ.1 ) CH11 = CHPOL
        IF ( MHT( IH ).EQ.2 ) CH11 = CHTOR
        IF ( MHT( IH ).EQ.3 ) CH11 = CHTHE
        IF ( MHC( IH ).EQ.1 ) CH3  = CHCOS
        IF ( MHC( IH ).EQ.2 ) CH3  = CHSIN
        IF ( IMODE.EQ.1 ) WRITE ( LU, 88 )
     1       CH11, MHL( IH ), MHM( IH ), CH3
        IF ( IMODE.EQ.1 ) WRITE ( LU, 89 )
        DO IRN = 1, NR
           RAD = RI + H*( IRN - 1 )
           IND = ( IRN - 1 )*NH + IH
           IF ( IMODE.EQ.1 ) WRITE ( LU, 90 ) RAD, RVEC1( IND ),
     1                      RVEC2( IND )
        ENDDO
 100    CONTINUE
      ENDDO
 88   FORMAT('#  ',A11,' L= ',I3,' M= ',I3,' ',A3)
 89   FORMAT('#  -----------------------------')
 90   FORMAT(f20.9,f20.9,f20.9)
C
      RETURN
      END
C*********************************************************************
