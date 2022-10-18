C*********************************************************************
C subroutine Single Harmonic Kinetic Energy Evaluation Routine *******
C            -      -        -       -      -          -       *******
C Steve Gibbons Thu Oct 28 08:51:22 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns 0.5d0 * \int_{volume} A.A dV where A is a single vector    C
C harmonic - the type being indicated by MHT( ih ).                  C
C                                                                    C
C This can be either magnetic energy or kinetic energy.              C
C                                                                    C
C The value for A = ( A_r, A_theta, A_phi ) is returned in DKE( 1 ). C
C                                                                    C
C The value for A = ( A_r, 0 , 0 ) is returned in DKE( 2 ).          C
C                                                                    C
C If MHT( ih ) does not correspond to a poloidal or toroidal         C
C vector harmonic, DKE( 1 ) and DKE( 2 ) are both returned zero.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IH        : Number of harmonic to be evaluated.                C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C MHT defines what each scalar function in a solution vector         C
C represents.                                                        C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DKE       : Dimension ( 2 ).                                   C
C                                                                    C
C                 DKE( 1 ) is returned with the whole contribution   C
C                 to the kinetic energy from harmonic IH             C
C                                                                    C
C                 DKE( 2 ) is returned with kinetic energy of the    C
C                 radial component of harmonic IH                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHKEER( IH, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1                   NDRVS, NDRVM, NFDCM, SV, XARR, DKE, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IH, NDCS, NR, INARR( * ), MHT( * ), MHL( * ), MHP( * ),
     1        NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION XARR( NR ), SV( * ), DKE( 2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, IHD, ITYPE, L, NR2, IS
      DOUBLE PRECISION COEF, DK, SQRLL1, PI, DERV( 2 ), D0F, D1F,
     1                 RAD, ER, ETOT, FAC
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR2 = INARR( 2 )
      IF ( NR2.NE.NR ) THEN
        PRINT *,' Subroutine SHKEER.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NR2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DKE( 1 ) = 0.0d0
      DKE( 2 ) = 0.0d0
C     .
      ITYPE = MHT( IH )
      L     = MHL( IH )
      IS    = MHP( IH )
      DK    = SQRLL1( L )
C     .
C     . First do poloidal harmonics
C     .
      IF ( ITYPE.EQ.1 .OR. ITYPE.EQ.4 ) THEN
        IHD = 1
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
          D1F = DERV( 2 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ER   = D0F*D0F
          ETOT = DK*DK*D0F*D0F + (D0F + RAD*D1F)*(D0F + RAD*D1F)
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
          DKE( 2 ) = DKE( 2 ) + ER*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
        DKE( 2 ) = DKE( 2 )*FAC*DK*DK
C       .
      ENDIF
C     .
C     . Now do toroidal harmonics
C     .
      IF ( ITYPE.EQ.2 .OR. ITYPE.EQ.5 ) THEN
        IHD = 0
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ETOT = RAD*RAD*D0F*D0F
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
