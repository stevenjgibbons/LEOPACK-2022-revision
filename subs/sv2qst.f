C*********************************************************************
C Subroutine Solution Vector 2 QST ***********************************
C            -        -      - --- ***********************************
C Steve Gibbons 17.6.99                                              C
C____________________________________________________________________C
C For a given radial node IRN;  SV2QST converts the poloidal and     C
C toroidal velocity coefficients from the SOLUTION VECTOR into QST   C
C coefficients.                                                      C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH	: Level of harmonics.                                C
C     NDIM      : Length of vector, SV.                              C
C     NR        : Number of radial grid nodes                        C
C     NH        : Number of harmonics (all types)                    C
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
C                                                                    C
C   ( but make a cautious note of input IPTF .... )                  C
C                                                                    C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     IVLMF     : Velocity / Magnetic field switch.                  C
C                                                                    C
C                  ivlmf = 1 --> put in velocity entries             C
C                  ivlmf = 2 --> put in magnetic field entries       C
C                                                                    C
C     IPTF      : The poloidal/toroidal flag.                        C
C                 This is set to 1 under the normal circumstances    C
C                 where the poloidal parts are found in MHT( IH )    C
C                 = 1 and the toroidal parts in MHT( IH ) = 2.       C
C                 However, if the curl is taken, then these are      C
C                 reversed and it is convenient to store them with   C
C                 the opposite definitions. In this case, set IPTF   C
C                 to 2.                                              C
C                                                                    C
C     IRN       : Number of radial node at which values to be        C
C                                                    taken.          C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SV        : Solution Vector Dimensions                         C
C                  Dim ( NDIM ) with NDIM = NH * NR                  C
C     RI        : Radius of inner core.                              C
C     RO        : Radius of outer core.                              C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector                           C
C                  Has dimensions (  LHMAX*(LHMAX+2) , 3).           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SV2QST ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC, IRN,
     1                    IVLMF, IPTF, ORD, SV, RI, RO, QST )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NDIM, NR, NH, MHT( NH ), MHL( NH ),
     1        MHM( NH ), MHC( NH ), IRN, IPTF, IVLMF
      DOUBLE PRECISION RI, RO, QST ( LH*(LH+2) , 3 ), SV( NDIM )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IOP, NDIM1, NDIM2, IH, L,
     1        M, ICS, INDSHC, NOHARM, IC, IPOL, ITOR
      DOUBLE PRECISION ZERO, POL, TOR, D0F, D1F, D2F, D3F, D4F,
     1                 DERIV, Q, SQRLL1, H, RAD, TOL
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check validity of IPTF
C
      IF ( IPTF.NE.1 .AND. IPTF.NE.2 ) THEN
        PRINT *,' Subroutine SV2QST, bad definition for IPTF '
        PRINT *,' IPTF = ',IPTF
        STOP
      ENDIF
C
      IF ( IVLMF.NE.1 .AND. IVLMF.NE.2 ) THEN
         PRINT *,' Subroutine SV2QST. IVLMF = ', IVLMF
         PRINT *,' Must be either one or two.'
         STOP
      ENDIF
C
      IF ( IVLMF.EQ.1 ) THEN
        IPOL = 1
        ITOR = 2
      ENDIF
C
      IF ( IVLMF.EQ.2 ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine SV2QST, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Set working parameters
      H = (RO-RI)/DBLE( NR - 1 )
      RAD = RI + H*DBLE( IRN - 1 )
C
C Check against RAD = 0; we divide by RAD later on ...
      IF (ABS(RAD).LT.TOL) THEN
         WRITE (6,987)
         WRITE (6,988)
         WRITE (6,989)
 987     FORMAT (' Subroutine SV2QST. The radius at JIRN ')
 988     FORMAT ('is too small and will cause a division ')
 989     FORMAT ('by zero error. Program aborted.')
         STOP
      ENDIF
C
C Check against too small a number of radial grid points
C
      IF ( NR .LT.8) THEN
         WRITE (6,887)  NR 
         WRITE (6,888)
         WRITE (6,889)
 887     FORMAT (' You are only using ',I2,' grid nodes. ')
 888     FORMAT (' I recommend using no less than 10; ')
 889     FORMAT (' Program aborted. ')
         STOP
      ENDIF
C
C Check against IRN > NR
      IF ( IRN.GT.NR ) THEN
         WRITE (6,789)
         WRITE (6,889)
 789     FORMAT (' IRN is greater than NR. ')
         STOP
      ENDIF
C ............ first let's set QST to zero ...........................
C
      NDIM1 = LH*(LH+2)
      NDIM2 = 3
      IOP = 0
      CALL MATOP ( QST, ZERO, NDIM1, NDIM2, IOP )
C
      DO IH = 1, NH
C       .
C       . First consider the case of a poloidal harmonic
C       .
        IF ( ( MHT( IH ).EQ.IPOL .AND. IPTF.EQ.1 ) .OR.
     1       ( MHT( IH ).EQ.ITOR .AND. IPTF.EQ.2 )       ) THEN
           L    = MHL( IH )
           M    = MHM( IH )
           ICS  = MHC( IH )
           NOHARM = INDSHC( L, M, ICS )
           Q = DBLE( L )
           CALL SVDERF ( NDIM, NH, IH, NR, IRN, SV, ORD,
     1                   D0F, D1F, D2F, D3F, D4F, RI, RO )
           POL = D0F
           DERIV = D1F + D0F/RAD
           QST( NOHARM, 1 ) = POL*Q*(Q+1.0d0)/RAD
           QST( NOHARM, 2 ) = DERIV*SQRLL1( L )
        ENDIF
C       .
C       . Now consider the case of a toroidal harmonic
C       .
        IF ( ( MHT( IH ).EQ.IPOL .AND. IPTF.EQ.2 ) .OR.
     1       ( MHT( IH ).EQ.ITOR .AND. IPTF.EQ.1 )       ) THEN
           L    = MHL( IH )
           M    = MHM( IH )
           ICS  = MHC( IH )
           NOHARM = INDSHC( L, M, ICS )
           IC = ( IRN - 1 )*NH + IH
           TOR = SV( IC )
           QST( NOHARM, 3) = -1.0d0*TOR*SQRLL1( L )
        ENDIF
C       .
      ENDDO
C
      RETURN
      END
C*********************************************************************

